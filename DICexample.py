#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is an example script that uses DIC data from Carrol et al
as input for SIF to find K field and cracktip data

The output is written to a CSV file 

@author: Swati Gupta
"""
import SIF_final as SIF
import numpy as np
from os import walk
import pdb
from datetime import datetime
import csv, timeit
from pandas import read_csv

synData = 0
r2 = 100
alpha = 2
N=-1
noise = -1

# provide path to directory containing displacement data files
path = 'DICresults/'
filenames = next(walk(path), (None, None, []))[2]
# exclude temp files 
for file in filenames:
    if file.startswith('.'):
        filenames.remove(file)
filenames.sort() #sort assuming files are named chronologically

#initilize
lenF = len(filenames)
cracktip1 = np.zeros((lenF-1,2)) # Geo
cracktip2 = np.zeros((lenF-1,2)) # Sep
cracktip3 = np.zeros((lenF-1,2)) # DC
K_field1 = np.zeros((lenF-1,6))  # Geo
K_field2 = np.zeros((lenF-1,6))  # Sepp
K_field3 = np.zeros((lenF-1,6))  # DC
discr = 10
geo = 0 # = 1 if use geometrical method too or 0 if only use separability

mat_constants = [0.327,109.1,43] #poisson's ratio, E, shear modulus
prev = [545, 315]
startT = timeit.default_timer()

#loop over the set of files
for i in range(0,lenF-1):
    file1 = path+filenames[i]
    print('filename \n', file1) 
    # pdb.set_trace()
    data = read_csv(file1)
    x = data['x'].values
    y = data['y'].values
    u = data['u'].values
    v = data['v'].values
    cracktip1[i],cracktip2[i],cracktip3[i],K_field1[i],K_field2[i], K_field3[i] = \
        SIF.SIF_projection(synData, r2 = 50, alpha = 2.5,coords = np.array([x,y]).T, 
                           coords_ref = np.array([x+u, y+v]).T, guess = prev, h=discr,geo = geo, constants = mat_constants)
    
    # prev = cracktip2[i]
    
endT = timeit.default_timer()
print("Time taken:", round(endT - startT, 2), "seconds to analyze ", lenF, "files")


## write output to file ##
currentDT = datetime.now() # get current date  and time
outputFile = "DIC_" + str(currentDT) + ".csv"
with open(outputFile, 'w') as f:
    writer = csv.writer(f)
    if geo:
        writer.writerow(['S.no','filename','x_geo','y_geo','K_I_geo','K_II_geo','T_geo','x_sep','y_sep','K_I_sep','K_II_sep','T_sep']) 
        writer.writerows(zip(range(1,lenF),filenames,cracktip1[:,0], cracktip1[:,1],K_field1[:,2],K_field1[:,3],K_field1[:,4],
                             cracktip2[:,0], cracktip2[:,1],K_field2[:,2],K_field2[:,3],K_field2[:,4]))
    else:   
        writer.writerow(['S.no','filename','x','y','K_I','K_II','T']) 
        writer.writerows(zip(range(1,lenF),filenames,cracktip2[:,0], cracktip2[:,1],K_field2[:,2],K_field2[:,3],K_field2[:,4]))
        
        
        
        