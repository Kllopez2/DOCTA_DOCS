# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 13:09:57 2025

@author: L093588
"""

import warnings
warnings.filterwarnings("ignore")
from scipy.spatial import ConvexHull, distance
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
import os
from glob import glob
import seaborn as sns


#--------------------Functions-------------------------
def select_files():
    #Function: pop-up window that allows you to select multiple files for analysis and saves the file paths. only opens TSV files but can be modified to open any type of file

    root = tk.Tk()
    root.withdraw()  # Hide the root window
    file_paths = filedialog.askopenfilenames(title="Select TSV files to analyze", filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")])

    return file_paths

def TSV_to_Slides(df):
    #Function: Finds the center point of each cell, drops unnecessary columns, finds the slide names, sorts all cells into respective slide dfs, stores dfs as a dictionary

    # Calculate the center of the two points for each x and two points for each y
    df['X'] = df[['Islet Object Info - Envelope left', 'Islet Object Info - Envelope right']].mean(axis=1)
    df['Y'] = df[['Islet Object Info - Envelope top', 'Islet Object Info - Envelope bottom']].mean(axis=1)

    # Drop the original columns
    df = df.drop(columns=['Islet Object Info - Envelope left', 'Islet Object Info - Envelope right', 'Islet Object Info - Envelope top', 'Islet Object Info - Envelope bottom'])

    # Get unique islets in the 'ROIType' column
    name_find = df['Name'].unique()

    # Create a dictionary to store DataFrames
    slides = {}

    # Loop through each unique islet and create a DataFrame
    for name in name_find:
        slides[f'{name}'] = df[df['Name'] == name]

    # print("Islet count:", np.max(islet_find))
    
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
    #Function: Removes islets with less than 4 glucagon/aldolase cells here it's named as "alpha", it's the name assigned w/ the AI software channel name

    #create a copy 
    fully_filtered_dfs =  filtered_dataframes.copy()
   
    for key, df in filtered_dataframes.items():

        #calculate the sum of A and B cells 
        sum_A_B = (df['Islet Object Info - LabelName'] == 'Beta').sum() + (df['Islet Object Info - LabelName'] == 'Alpha').sum()
        if sum_A_B < 10:
            del fully_filtered_dfs[key]
            
    return fully_filtered_dfs

   
def cell_types_coordinates(working_islet):
    #Function: Separates the cell types into individual dfs, creates a convex hull of the points, and finds the area of the convex hull

    #beta and alpha pure cell coordinates
    aldo_coordinates = working_islet[working_islet['Islet Object Info - LabelType'] == 6][['X', 'Y']].values
    beta_coordinates = working_islet[working_islet['Islet Object Info - LabelType'] == 5][['X','Y']].values
    other_cells = working_islet[working_islet['Islet Object Info - LabelType'] == 2][['X','Y']].values

    #adding other cells to calculate the convex hull
    beta_aldo_coordinates = np.vstack((aldo_coordinates, beta_coordinates, other_cells))
    
    aldo_cells_ratio = len(aldo_coordinates) / (len(beta_coordinates) +len(aldo_coordinates))

    #Calculating the convex hull using the x,y coordinates of the beta cells
    islet_hull = ConvexHull(beta_aldo_coordinates)
    area = islet_hull.volume *1000000
    
    
    
    return aldo_cells_ratio, area


def trim_means(all_data):
    #Function: Stores only the values between the low and high % specificications 

    # Calculate Q1 (25th percentile) and Q3 (75th percentile)
    #filtered_df = islet_analysis_image_results[(islet_analysis_image_results['Median Length/Area Value - 1/μm'] != 0)]
    lower_percentile = all_data['Islet Area - μm²'].quantile(0.05)
    upper_percentile = all_data['Islet Area - μm²'].quantile(0.95)
    df_no_outliers = all_data[(all_data['Islet Area - μm²'] >= lower_percentile) & (all_data['Islet Area - μm²'] <= upper_percentile)]
    
    return df_no_outliers


def plot_distribution_histogram(alpha_cells_histogram, measurement, xlim, binning):
    #Functions: Normal histogram plotter but made easier to include within the script

    plt.hist(alpha_cells_histogram[measurement], bins= binning, edgecolor='black')
    plt.xlim(0, xlim)
    plt.xlabel(measurement)
    plt.ylabel('Counts')
    plt.title(experimental_conditions)
    plt.show()
    
    return 



def beta_cells_density_graphs(beta_cells_ratio):
    #Function: Creates scatter plots for the islet area vs beta cell ratio information 
        
    plt.scatter(beta_cells_ratio['Islet Area - μm²'], beta_cells_ratio['Beta Cells Ratio'], color = 'green')
    plt.title(experimental_conditions)
    plt.xlabel('Islet Area (μm²)')
    plt.ylabel('Beta Cells Density')
    
    plt.ylim(0, 1)
    plt.xlim(0, 150000)
    plt.show()
    

#------------------------------------main script-------------------------------



file_paths = select_files()

all_data = pd.DataFrame()
df = pd.DataFrame() 
final_islet_count = []   
alpha_big_histogram = pd.DataFrame()

for file_path in file_paths:
    
    glucose_value = int(os.path.basename(file_path[-7:-4]))

    #----------------------Open the File----------------- previous script
    
    # Read the TSV file into a DataFrame
    df = pd.read_csv(file_path, sep='\t') #here is where you should collect the file name later on!!
    # Grab a single value
    day = df.loc[1, 'Study level 3']  # This grabs the day 
    condition = df.loc[1, 'Study level 4'] #This grabs the condition (HFD vs Control)
    animal = df.loc[1, 'Name'] #grabbing slide name 
    experimental_conditions = condition + ' ' + day 
    
    print(f"----------------Running File {animal}------------------") #for my reference!
    print(f"Currently Running: {experimental_conditions}")
        
    # Drop the first 5 columns
    df = df.drop(df.columns[:4], axis=1) 
               
    slides = TSV_to_Slides(df) #now this returns all slides separated by name. For each name, Separate by Islets
    islets = slides_to_islets(slides)


    
    initial_islet_count = len(islets)
            
    fully_filtered_dfs = remove_islets(islets)    
            
    islet_analysis_image_results = pd.DataFrame()
    
    filtered_islet_count = len(fully_filtered_dfs)
    lost_islets = initial_islet_count - filtered_islet_count
    print(f"Final Islet Count: {filtered_islet_count} ")
    print(f"Islets Lost: {lost_islets} ")
    
    #going through each individual islet here to extract data
    for label, df in fully_filtered_dfs.items():
      
        name = df.iloc[0, df.columns.get_loc('Name')]

        aldo_cells_ratio, area = cell_types_coordinates(df)
        
        temp_df = pd.DataFrame({
            'Animal': [name],
            'Islet Label':[label],
            'Islet Area - μm²': [area],
            'Glucose' : [glucose_value],
            'AldoB Cells Ratio': [aldo_cells_ratio], 
            })
        
        #adding the current analyzed dataframe to the large analysis file
        islet_analysis_image_results = pd.concat([islet_analysis_image_results, temp_df], ignore_index = True)
        

    #df_no_outliers = islet_analysis_image_results.dropna(subset=['Median Length/Area Value - 1/μm'])
    final_islet_count.append(filtered_islet_count) 
    
    all_data = pd.concat([all_data, islet_analysis_image_results], ignore_index = True)


all_data_no_outliers = trim_means(all_data)

total_islets = sum(final_islet_count)

plot_distribution_histogram(all_data,'Islet Area - μm²', 18000000, 50)
plot_distribution_histogram(all_data,'AldoB Cells Ratio', 1, 25)

plot_distribution_histogram(all_data_no_outliers,'Islet Area - μm²', 200000, 50)
plot_distribution_histogram(all_data_no_outliers,'AldoB Cells Ratio', 1, 25)

#Get the directory of the first selected file to save the output in the same folder------ real output files
file_name2 = experimental_conditions.replace(" ", "_") + ".xlsx"
output_folder = os.path.dirname(file_path)
all_dataoutput_file = os.path.join(output_folder, file_name2)
all_data_no_outliers.to_excel(all_dataoutput_file, index=False)

