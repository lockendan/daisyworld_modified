# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 12:51:49 2020

@author: Daniel
"""
import numpy as np


# Define constants and variables
alphaw = 0.01   # Cover fraction of white daisies
alphab = 0.01   # Cover fraction of black daisies
p = 1           # The fraction of habitable surface
alphag = p-alphaw-alphab # Cover fraction of bare ground
aw = 0.75       # Albedo of white daisies
ab = 0.25       # Albedo of black daisies
ag = 0.5        # Albedo of bare ground
gamma = 0.3     # The death rate 
S = 1000        # Solar constant (W/m^2)
maxn = 1000     # Maximum number of iterations
tol = 0.000001  # Tollerance of solution
luminosities = np.arange(0.5,1.6, 0.002) # Stelar luminosities
alphaw_out = np.ones(len(luminosities))*np.nan # Output variable for white
alphab_out = np.ones(len(luminosities))*np.nan # Output variable for black
temp_out = np.ones(len(luminosities))*np.nan   # Output variable for temp
