# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 13:09:57 2025

@author: Karen Leonor Lopez - 2025 Spring TDA Intern
"""

import warnings
warnings.filterwarnings("ignore")

import openpyxl
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull, distance, KDTree
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import os
from glob import glob
import seaborn as sns
from scipy.spatial import Delaunay, KDTree


"""
------------------------------functions start here-------------------------------
Note: All functions are in the order in which they're used!
"""

def select_files():
    #Function: pop-up window that allows you to select multiple files for analysis and saves the file paths. only opens TSV files but can be modified to open any type of file
    
    root = tk.Tk()
    root.withdraw() 
    file_paths = filedialog.askopenfilenames(title="Select TSV files to analyze", filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")])
       
    return file_paths 


def TSV_to_Slides(df):
    #Function: Finds the center point of each cell, drops unnecessary columns, finds the slide names, sorts all cells into respective slide dfs, stores dfs as a dictionary
    
    df['X'] = df[['Islet Object Info - Envelope left', 'Islet Object Info - Envelope right']].mean(axis=1)    # Calculate the center of the two points for each x and two points for each y
    df['Y'] = df[['Islet Object Info - Envelope top', 'Islet Object Info - Envelope bottom']].mean(axis=1)
    df = df.drop(columns=['Islet Object Info - Envelope left', 'Islet Object Info - Envelope right', 'Islet Object Info - Envelope top', 'Islet Object Info - Envelope bottom'])    # Drop the original columns
    name_find = df['Name'].unique()    # Get unique islets in the 'ROIType' column
    slides = {}    # Create a dictionary to store DataFrames
    for name in name_find:    # Loop through each unique islet and create a DataFrame
        slides[f'{name}'] = df[df['Name'] == name]
    
    return slides


def slides_to_islets(slides):
   #Function: Assigns a new number to all islets to prevent having islets from different slices w/ same numbers. Saves islets as individual dfs. 
    islets = {}
 
    counter = 0
    for key, df in slides.items():
        islet_find = df['Islet Object Info - ROIType'].unique()
        
        for islet in islet_find:
             counter += 1
             islets[f'{counter}'] = df[df['Islet Object Info - ROIType'] == islet]
                                   
    return islets


def remove_islets(filtered_dataframes):
    #Function: Removes islets with less than 4 Alpha cells (so any islet without Alpha cells too)
    #create a copy 
    fully_filtered_dfs =  filtered_dataframes.copy()
   
    for key, df in filtered_dataframes.items():

        sum_Alpha = (df['Islet Object Info - LabelName'] == 'Alpha').sum() 
        if sum_Alpha < 4:
            del fully_filtered_dfs[key]
            
    return fully_filtered_dfs


def cell_types_coordinates(working_islet):
    #Function: Separates the cell types into individual dfs, creates a convex hull of the points, and finds the area of the convex hull
    
    alpha_coordinates = working_islet[working_islet['Islet Object Info - LabelType'] == 6][['X', 'Y']].values
    beta_coordinates = working_islet[working_islet['Islet Object Info - LabelType'] == 5][['X','Y']].values
    other_cells = working_islet[working_islet['Islet Object Info - LabelType'] == 2][['X','Y']].values

    beta_alpha_coordinates = np.vstack((alpha_coordinates, beta_coordinates, other_cells))    #adding all cells to calculate the convex hull
    
    other_cells_ratio = len(other_cells) / (len(other_cells)+ len(beta_coordinates) +len(alpha_coordinates))     #finds the ratio of whichever cell of interest 
    alpha_cell_count = len(alpha_coordinates)    #counts number of alpha cells in the islet
    islet_hull = ConvexHull(beta_alpha_coordinates)     #Calculates the convex hull using the function below 
    
    return alpha_coordinates, islet_hull, beta_alpha_coordinates, beta_coordinates, other_cells_ratio, alpha_cell_count, other_cells

def distance_to_hull(alpha_coordinates, islet_hull):
    #Function: Finds the distance from each point to the nearest part of the convex hull, saves the distance

    distances_to_hull = []
    for point in alpha_coordinates:
        min_distance = float('inf')
        nearest_point = None
        for simplex in islet_hull.simplices:
            segment_start = beta_alpha_coordinates[simplex[0]]
            segment_end = beta_alpha_coordinates[simplex[1]]
            distance, nearest = point_to_segment_distance(point, segment_start, segment_end)
            if distance < min_distance:
                min_distance = distance
                nearest_point = nearest
        distances_to_hull.append((min_distance, nearest_point))

    distance_df = pd.DataFrame(distances_to_hull, columns = ['Alpha Cell Distance', 'Nearest Point'])     #Saves distances as a df
    
    return distance_df

def point_to_segment_distance(point, segment_start, segment_end):
    #Function: Finds the nearest point on the convex hull to each alpha cell, creates a line when you plot it from the cell to the convex hull 
        
    line_vec = segment_end - segment_start
    point_vec = point - segment_start
    line_len = np.linalg.norm(line_vec)
    line_unitvec = line_vec / line_len
    point_vec_scaled = point_vec / line_len
    t = np.dot(line_unitvec, point_vec_scaled)
    t = np.clip(t, 0, 1)
    nearest = segment_start + t * line_vec
    distance = np.linalg.norm(point - nearest)
    return distance, nearest


def cell_displacement(islet_hull, alpha_coordinates):
    #Function: Creates a df with the name, islet label, cell to convex hull centroid, distance from cell to convex hull, and islet area
    
    centroid = np.nanmean(beta_alpha_coordinates[islet_hull.vertices], axis = 0, keepdims=True) #center point of the convex hull
    
    area = islet_hull.volume *1000000   #Finds the convex hull area 

    alpha_displacement = pd.DataFrame()
    
    for (alpha_coordinate, distance1) in zip(alpha_coordinates, distance_df['Alpha Cell Distance']): #Loops used to calculate the distance from each point to the centroid
        
        D1 = np.linalg.norm(alpha_coordinate - centroid) *1000         #calculate the distance from the centroid to alpha coordinate
            
        D2 = distance1 * 1000        #distance from alpha coordinate to the convex hull

        new_row = pd.DataFrame({
            'Animal': [animal],
            'Islet Label': [label],
            'Cell to Centroid': [D1], 
            'Cell to ConvexHull': [D2], 
            'Islet Area - μm²': [area],
            })
        
        alpha_displacement = pd.concat([alpha_displacement, new_row], ignore_index=True)
        
    return alpha_displacement, area


def trim_means(islet_analysis_image_results, measurement, low, high):
    #Function: Stores only the values between the low and high % specificications 
    
    lower_percentile = islet_analysis_image_results[measurement].quantile(low)
    upper_percentile = islet_analysis_image_results[measurement].quantile(high)
    df_no_outliers = islet_analysis_image_results[(islet_analysis_image_results[measurement] >= lower_percentile) & (islet_analysis_image_results[measurement] <= upper_percentile)]
    
    return df_no_outliers


def plot_distribution_histogram(alpha_cells_histogram, measurement, xlim, binning):
    #Functions: Normal histogram plotter but made easier to include within the script
    
    plt.xlim(0, xlim)
    plt.hist(alpha_cells_histogram[measurement], bins= binning, edgecolor='black')
    plt.xlabel(measurement)
    plt.ylabel('Frequency')
    plt.title(experimental_conditions)
    plt.show()
    return 



"""
-------------------------The Main Script Starts Here------------------------------
"""


file_paths = select_files() #Opens a pop-up window to select the files you want to analyze

#making the dataframes I want to use in the for loop
all_data = pd.DataFrame()
all_alpha_cells_histogram = pd.DataFrame()
alpha_cells_histogram = pd.DataFrame() 
final_islet_count = []   
#all_cells_count = []

for file_path in file_paths:
    #Function: loops though each file path 
    glucose_value = int(os.path.basename(file_path[-7:-4])) #grabs the glucose value, which i save next to the animal name. Example: A0231_120 <--- last 3 digits
    
    df = pd.read_csv(file_path, sep='\t') #using the first file path, opens the file, saves it as a dataframe
    day = df.loc[1, 'Study level 3']  # This grabs the day saved inside the original TSV structure 
    condition = df.loc[1, 'Study level 4'] #This grabs the condition (HFD vs Control) saved inside the original TSV structure
    animal = df.loc[1, 'Name'] #grabbing slide name saved inside the original TSV structure
    experimental_conditions = condition + ' ' + day 
    
    print(f"----------------Running File {animal}------------------") #for reference, to verify all files are run 
    print(f"Currently Running: {experimental_conditions}") 
        
    df = df.drop(df.columns[:4], axis=1) #Drops the first four columns, they are not necessary from here on
               
    slides = TSV_to_Slides(df) #For each animal file, converts each slice into its own dataframe, stores all cell values into their respective slices, and goes from 4 points to a single point. 
    
    islets = slides_to_islets(slides) #Loops through all slides and assigns a new number to each islet to prevent same-number islets
    
    initial_islet_count = len(islets) #counts initial # of islets to see how many are lost w/ filters
            
    fully_filtered_dfs = remove_islets(islets) #Removes islets with less than 4 alpha cells
    
    islet_analysis_image_results = pd.DataFrame()
    
    filtered_islet_count = len(fully_filtered_dfs) 
    lost_islets = initial_islet_count - filtered_islet_count 
    print(f"Final Islet Count: {filtered_islet_count} ") #For reference that code is running properly 
    print(f"Islets Lost: {lost_islets} ")
    
    for label, df in fully_filtered_dfs.items(): 
        #Function: Loops through each islet df to calculate alpha cell distances 
        name = df.iloc[0, df.columns.get_loc('Name')] #Grabs the slide name 

        alpha_coordinates, islet_hull, beta_alpha_coordinates, beta_coordinates, other_cells_ratio, alpha_cell_count, other_cells = cell_types_coordinates(df) #for each islet, separates the alpha, beta, and other cells; counts # of alpha cells; calculates the convex hull
        
        distance_df = distance_to_hull(alpha_coordinates, islet_hull) #finds the distances from the alpha cells to the convex hull
        
        results, area = cell_displacement(islet_hull, alpha_coordinates) #Grabs the convex hull area, sorts parameters 
        
        #all_cells_count.append(len(beta_alpha_coordinates))
        
        index_value = np.nanmedian(results['Cell to ConvexHull']) #For each islet, finds the median distance from convex hull to alpha cell
        
        temp_df = pd.DataFrame({
            'Animal': [name],
            'Glucose Value': [glucose_value],
            'Islet Label': [label],
            'Islet Area - μm²': [area],
            'Alpha Cell Count': [alpha_cell_count],
            'Other Cells Ratio': [other_cells_ratio],
            'Islet Fragmentation Index': [index_value]
            })
        
        islet_analysis_image_results = pd.concat([islet_analysis_image_results, temp_df], ignore_index = True) #saves single set of results from each islet
        
        alpha_cells_histogram = pd.concat([alpha_cells_histogram, results], ignore_index = True)    #saves ALL alpha cell distances     

    all_data_no_outliers = islet_analysis_image_results.dropna(how = 'any')     #removes any islet with any zero (there really shouldn't be any), saved results for all islets in one file
    all_alpha_cells_histogram = pd.concat([all_alpha_cells_histogram, alpha_cells_histogram], ignore_index = True)      #saves all alpha cells within the one file
    final_islet_count.append(filtered_islet_count) #Counts the number of islets analyzed per condition
    
    all_data = pd.concat([all_data,  all_data_no_outliers], ignore_index = True) #Saves all values for all animals 
    
    
"""
--------------------------------Loops end here-----------------------------------
"""
print(f"Total Islets Analized for {experimental_conditions}: {sum(final_islet_count)} ")
print(f"Total Alpha Cells Analyzed: {len(all_alpha_cells_histogram['Cell to ConvexHull'])}")

trimed_alpha_cells = trim_means(all_alpha_cells_histogram, 'Islet Area - μm²', 0.02, 0.98) #trims the top and lowest 2% to remove outliers based on area - my concern is overly large islets made from multiple islets
trimed_all_data = trim_means(all_data, 'Islet Area - μm²', 0.02, 0.98) #repeated for the median values

"""
----------------------------Saving/Exporting the Files-------------------------------
"""

#Get the directory of the first selected file to save the output in the same folder
file_name = experimental_conditions.replace(" ", "_") + "all_islets" + ".xlsx" #creates export file name for all islet medians 
file_name2 = experimental_conditions.replace(" ", "__") + "all_cells" + ".xlsx" #creates export file name for all alpha cell values
output_folder = os.path.dirname(file_path) 
output_file = os.path.join(output_folder, file_name)
output_file2 = os.path.join(output_folder, file_name2)
trimed_all_data.to_excel(output_file, index=False) #exporting all islets to excel
trimed_alpha_cells.to_excel(output_file2, index=False) #exporting all alpha cell values to excel 

plot_distribution_histogram(trimed_alpha_cells,'Cell to ConvexHull', 200, 50) #plots the histograms to let me know it's all done!
plot_distribution_histogram(trimed_all_data,'Islet Area - μm²', 30000000, 300)
plot_distribution_histogram(trimed_all_data,'Islet Fragmentation Index', 200, 100)


#Exporting the sheets

# Slide1 = slides['R0301 09']

# Slide1.to_excel(output_file2, index = False)
