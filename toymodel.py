#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 08:16:28 2023

@author: swati gupta

ls

This is an example script that creates synthetic data for 
SIF to find K field and cracktip data 
"""
import SIF_final as SIF
import numpy as np

choice = 1
r2 = 100
alpha =4
N= 50000
noise = 0
trials = 1
discr = 2
geo = 1
shape = 3
Dc = 1

cracktip1 = np.zeros((trials,2)) # Geometrical (geo) method
cracktip2 = np.zeros((trials,2)) # Separability (Sep)
cracktip3 = np.zeros((trials,2)) # Displacement Correlation (DC)
K_field1 = np.zeros((trials,6))  # Geo
K_field2 = np.zeros((trials,6))  # Sep
K_field3 = np.zeros((trials,6))  # Dc
            

for i in range(0,trials):
    cracktip1[i],cracktip2[i],cracktip3[i],K_field1[i], K_field2[i], K_field3[i] = \
    SIF.SIF_projection(choice,shape,N, r2, alpha,noise,guess = [5,5], h=discr,geo = 0, PS=0, DC = 0)
 