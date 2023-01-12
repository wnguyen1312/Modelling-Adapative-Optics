#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This file is used to optimise a bi-convex lens set up for a given object and
image distance.

Tips: 
    - Use optimum.x to call the optimised parameters 
    - Use optimum.fun to call the optimised RMS spot size value 

Created on Mon Nov  7 10:49:21 2022

@author: williamnguyen
"""

import optimisation 
import scipy.optimize as op 


c0 = 0.02 # initial guess for the first curvature 
c1 = -0.02 # initial guess for second curvature 

# sets focus to the paraxial focus of the planoconvex lens to compare rms.
focus_point = 198.452700178 


#setting up bound for optimisation

rmin = 5 #has to be at least 5 cm in radius for 10 mm beam 


#can define rmax to impose max limit on lens size if desired. 


bound = ((0,1./(rmin)),(-1./(rmin),0)) #one convex side, one concave side


optimum = op.minimize(optimisation.system_rms, [c0, c1], args = (100, 105,
                                     focus_point, 10 ), bounds=bound)

optimum_curvature = optimum.x
min_rms = optimum.fun

print('The optimum curvature combination is:', optimum_curvature
      , 'mm^-1')
print('The minimised RMS is:', min_rms, 'mm')
        
#Performance is better than planoconvex for 10 mm bundle 

