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

def find_medians(library):
    
    library_medians = pd.DataFrame()
    
    for name, df in library.items():
        
        area = np.median(df['Islet Area - μm²'])
        fragmentation_index = np.median(df['Islet Fragmentation Index'])
        name2 = df.iloc[0, df.columns.get_loc('Animal_5')]
        glucose_value = df.iloc[0, df.columns.get_loc('Glucose Value')]
        median_ratio = np.median(df['Other Cells Ratio'])

        new_row = pd.DataFrame({
             'Animal': [name2],
             'Glucose' : [glucose_value],
             'Fragmentation Index': [fragmentation_index],
             'Area': [area],
             'Other Cells Ratio': [median_ratio],
             })
         
        library_medians = pd.concat([library_medians, new_row], ignore_index=True)        
        
    return library_medians

def average_means(library_medians):
    
    average_medians = pd.DataFrame()
    
    # Iterate through each unique name
    grouped = library_medians.groupby('Animal')    
    
    averages_df = grouped.mean().reset_index()
    
    return averages_df





#---------------------code starts here-------------------------

file_paths = select_files()
df = pd.read_excel(file_paths)
condition = file_paths[118:125]


log_df = df

#log_df['Log Index'] = np.log(df['Cell to ConvexHull'])

library = slice_sorter(log_df)

library_medians = find_medians(library)

library_averages = average_means(library_medians)

#Get the directory of the first selected file to save the output in the same folder
file_name = condition + "_Animal_Averages_Islets" + ".xlsx"


output_folder = os.path.dirname(file_paths)


output_file2 = os.path.join(output_folder, file_name)
library_averages.to_excel(output_file2, index=False) #exporting all alpha cell values

