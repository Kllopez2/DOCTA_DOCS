# -*- coding: utf-8 -*-
"""
Created on Thu Mar  6 15:34:31 2025

@author: L093588
"""
import warnings
warnings.filterwarnings("ignore")

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
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks


def select_files():
    root = tk.Tk()
    root.withdraw()  # Hide the root window
    file_paths = filedialog.askopenfilename(filetypes=[("Excel files", "*.xlsx *.xls")])
    
    return file_paths

def slice_sorter(df):
    # Get unique Animals in the 'Animal' column
    df['Animal_5'] = df['Animal'].apply(lambda x: x[:5])

    animal_find = df['Animal'].unique()
    # Create a dictionary to store DataFrames
    library = {}
    # Loop through each unique slice and create a DataFrame
    for animal in animal_find:
        library[f'{animal}'] = df[df['Animal'] == animal]
    
    return library

def find_means(library):
    
    library_means = pd.DataFrame()
    
    for name, df in library.items():
        
        mean = np.mean(df['AldoB Cells Ratio'])
        area = np.mean(df['Islet Area - μm²'])
        STD = np.std(df['AldoB Cells Ratio'])
        glucose_value = df.iloc[0, df.columns.get_loc('Glucose')]
        name2 = df.iloc[0, df.columns.get_loc('Animal_5')]
        
        new_row = pd.DataFrame({
             'Animal': [name2],
             'Glucose Value': [glucose_value],
             'STD of Ratio': [STD], 
             'Area': [area],
             'AldoB Ratio': [mean],
             })
         
        library_means = pd.concat([library_means, new_row], ignore_index=True)        
        
    return library_means

def median_of_means(library_means):
    
    median_means = pd.DataFrame()
    
    #library_means['Animal'] = library_means['Animal'].str[:-2]
    # Iterate through each unique name
    grouped = library_means.groupby('Animal')    
    
    median_means = grouped.mean().reset_index()
    
    return median_means

def plot_distribution_histogram(alpha_cells_histogram, measurement, xlim, binning):
    
    plt.hist(alpha_cells_histogram[measurement], bins= binning, edgecolor='black')
    plt.xlim(0, xlim)
    plt.xlabel(measurement)
    plt.ylabel('Counts')
    plt.title(condition)
    plt.show()
    
    return 

# def gaussean_kde(df, measurement):
    
#     ax = sns.kdeplot(df['AldoB Cells Ratio'], bw_adjust=0.5)
#     x = ax.lines[0].get_xdata()
#     y = ax.lines[0].get_ydata()
    
#     peaks, _ = find_peaks(y, height = 0)
#     print(f'number of peaks: {len(peaks)}')
    
#     plt.plot(x, y)
#     plt.plot(x[peaks], y[peaks], "x")
#     plt.xlabel('AldoB Ratio')
#     plt.ylabel('Density')
#     plt.title(condition +'KDE Plot with Peaks')
#     plt.show()
    
#     return 


#---------------------code starts here-------------------------

file_paths = select_files()
df = pd.read_excel(file_paths)
condition = file_paths[132:142]


#log_df = df
#gaussean_kde(df, 'AldoB Cells Ratio')
#peaks = find_KDE_peaks()

library = slice_sorter(df)

library_means = find_means(library)

animal_medians = median_of_means(library_means)

#Get the directory of the first selected file to save the output in the same folder
file_name = condition + "_Animal_Averages_allCells" + ".xlsx"

plot_distribution_histogram(df, 'AldoB Cells Ratio', 1, 25)
plot_distribution_histogram(library_means,'AldoB Ratio', 1, 25)
plot_distribution_histogram(animal_medians,'AldoB Ratio', 1, 25)

output_folder = os.path.dirname(file_paths)

file_name = condition + "_Animal_Medians" + ".xlsx"
output_file1 = os.path.join(output_folder, file_name)
animal_medians.to_excel(output_file1, index=False) #exporting all alpha cell values

file_name2 = condition + "_Animal_slideMeans" + ".xlsx"
output_file2 = os.path.join(output_folder, file_name2)
library_means.to_excel(output_file2, index=False) #exporting all alpha cell values

