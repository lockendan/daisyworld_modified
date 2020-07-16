# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 12:53:33 2020

@author: Daniel
"""

import matplotlib.pyplot as plt

from main import *
# Plot the results
# Cover fractions
white = plt.plot(luminosities,alphaw_out*100,'b', label='White')
black = plt.plot(luminosities,alphab_out*100,'k', label='Black')
plt.legend(loc='upper right')
plt.xlabel('Luminosity')
plt.ylabel('Surface cover %')
plt.title('Cover fractions')
plt.show()

# Planetary temperature
plt.figure()
plt.plot(luminosities,temp_out-273.15,'r')
plt.xlabel('Luminosity')
plt.ylabel('Temperature (C)')
plt.title('Planetary temperature')
plt.show()