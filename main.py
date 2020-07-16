#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 19:19:16 2019

@author: tammas loughran
@email: t.loughran@lmu.de

This is an implementation of daisy world: a theoretical 
model of the Gaia hypothesis. This theory was initially developed by 
Lovelock (1983) to demonstrate the plausibility of living things 
interacting with, and regulating, their environment.

    Wood, A. J., G. J. Ackland, J. G. Dyke, H. T. P. Williams, and T. M. 
        Lenton, 2015: Daisywolrd: a Review. Rev. Geophys., 46, RG1001, 
        https://doi.org/10.1029/2006RG000217.

    Copyright (C) 2019  Tammas Loughran

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""


import numpy as np


from functions import *
from variables import *


if __name__=='__main__':

    # Main loop for changing luminosity
    for i,L in enumerate(luminosities):
        # Set a minimum for cover fractions
        if alphaw<0.01: alphaw = 0.01
        if alphab<0.01: alphab = 0.01
        alphag = p-alphaw-alphab
        # Reset counters
        n = 0
        changew, changeb = 1,1
        # Run loop for daisy earth.
        while (n<maxn) and (changew>tol) and (changeb>tol):
            # Store the initial cover fractions
            sw,sb = alphaw, alphab
            # Planetary albedo
            planet_albedo = albedo(alphaw,alphab,alphag,aw,ab,ag)
            # Planetary temperature
            T = planetary_temp(S,planet_albedo, L=L)
            # Local temperature
            Tw = local_temp(planet_albedo,aw,T)
            Tb = local_temp(planet_albedo,ab,T)
            # Birth rate
            betaw = beta(Tw)
            betab = beta(Tb)
            # Change in daisies
            dawdt = daisy_replicator(alphaw, alphag, betaw, gamma)
            dabdt = daisy_replicator(alphab, alphag, betab, gamma)
            # Integrate
            alphaw = euler(alphaw, dawdt)
            alphab = euler(alphab, dabdt)
            alphag = p-alphaw-alphab
            n += 1
        # Store the output
        alphaw_out[i] = alphaw
        alphab_out[i] = alphab
        temp_out[i] = T
    
    # Plot the results
    from plot_results import *
