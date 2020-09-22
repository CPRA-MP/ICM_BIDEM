# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 16:46:35 2020

@author: ddileonardo

Interpolate the output from BI-DEM, which has grid points that move, to a fixed 
grid for the ICM
"""

import pandas as pd #used here for reading csv files
import numpy as np #used for creating arrays of the csv files
from scipy.interpolate import griddata #used for the interpolation

#Define the save location
save_loc = 'xyz output name' #output name and file path

#Read the data files
model_results = 'model results filename' #This is the filepath and filename of the model results
results = pd.read_csv(model_results,delim_whitespace=True,names=['x','y','z','trans']) #Read the model results file

fixed_output = 'fixed xy filename' #This is the filepath and filename of the fixed output grid points
fixed_grid = pd.read_csv(fixed_output) #Read the csv file of the fixed output points; 1 headerline; format: X-coordinate, Y-coordinate, transect, point number

#Create arrays for griddata
results_xy = np.vstack((results.x,results.y)).T #array of xy model results (model grid)
results_z = np.array(results.z) #array of model z values 

fixed_grid_array = np.vstack((fixed_grid.XCoord,fixed_grid.YCoord)).T  #array of fixed output xy locations

#Interpolation 
interp_results = griddata(results_xy,results.z,fixed_grid_array,method='linear') #interp_results is the z value of the model results interpolated to the fixed output grid point

#Write the x, y, z values from the fixed grid
#File format is: point index, X, Y, Z value
#There is one headerline
#There is one file for each BI-DEM Region
write = pd.DataFrame(data = np.hstack((fixed_grid_array,np.expand_dims(interp_results,axis=1))),columns=['x','y','z'])
write.to_csv(save_loc)
