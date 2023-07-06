#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 06 12:46:21 2022
Modified on Mon, Jul 3, 2023

@author: Swati Gupta
"""

import SIF_final as SIF
import time
import numpy as np
from os import walk
import matplotlib.pyplot as plt
import read_lmp_file


# INPUT: path to dump files

def calcSIFs(mypath):  
    
    global coords_ref
    
    filenames = next(walk(mypath), (None, None, []))[2]
    filenames.sort()
    
    lenF = len(filenames)
    cracktip1 = np.zeros((lenF-1,2)) # Geo
    cracktip2 = np.zeros((lenF-1,2)) # Sep
    cracktip3 = np.zeros((lenF-1,2)) # DC
    K_field1 = np.zeros((lenF-1,6))  # Geo
    K_field2 = np.zeros((lenF-1,6))  # Sep
    K_field3 = np.zeros((lenF-1,6))  # DC
    
    pdb.set_trace()
    
    prev = [520,0] #starting crack tip location guess
    lmpRef = 0 #dummy variable to indicate if lmp ref object is being passed
    discr = 2 #discretization size for Separability
    confg = 0 # config = 1 when using file = config_pacman
    fileType = 'dump'
    const = [0.2832,71.8,25.81]
    
    startT = time.time()
    for i in range(0,lenF-1):
        
        file1 = mypath+filenames[i+1]
        file2 = mypath+filenames[i]
        
        coords, coords_ref = getCoords(file1,file2,fileType,lmpRef, confg)
        
        cracktip1[i],cracktip2[i],cracktip3[i],K_field1[i],K_field2[i], K_field3[i] = \
        SIF.SIF_projection(synData = 2, r2 = 60, alpha = 2.5, coords = coords, 
                           coords_ref = coords_ref,guess = prev, h=discr,
                           geo = 0,Sep = 1, DC = 0, constants = const)
        # prev = cracktip1[i]
        prev = [prev[0]- (0.2 * i), prev[1]]
        lmpRef = 1
        # confg = 1
        
    
    endT = time.time()
    print("Time taken:", endT-startT, "seconds")
        
    return cracktip1,cracktip2,cracktip3,K_field1,K_field2, K_field3

def getCoords(file1,file2,fileType,lmp_ref, config):
    global coords_ref
    
    lmp = read_lmp_file.read_lmp_file(filename = file1, filetype = fileType).read_lmp_file()
    if lmp_ref == 0:
        ref = read_lmp_file.read_lmp_file(filename = file2, filetype = fileType).read_lmp_file()
        
    if config:
        coords = lmp.coordinate
        coords_ref = lmp_ref.coordinate
    else:
        ids_SiO = np.where((lmp.type == 1) | (lmp.type == 2))[0]
        coords = np.array([lmp.coordinate[ids_SiO,1],lmp.coordinate[ids_SiO,0]]).T
        if lmp_ref == 0:
            ids_SiO_ref = np.where((ref.type == 1) | (ref.type == 2))[0]
            coords_ref = np.array([ref.coordinate[ids_SiO_ref,1],ref.coordinate[ids_SiO_ref,0]]).T
          
      
    return coords, coords_ref