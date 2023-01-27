#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is used to model adaptive optics by using aspheric surfaces: 

1. Optimised single aspheric surface (compare with single spherical surface)

2. Optimised plano-aspheric lens (compare with plano-convex lens)

3. Optimised bi-convex lens by using aspheric surface (compare with bi-convex)

Tips: 
- Use optimum.x to call the optimised parameters 
- Use optimum.fun to call the optimised RMS spot size value 

Created on Sun Nov 20 13:13:25 2022

@author: williamnguyen
"""


import scipy.optimize as op
import optimisation 

#%% 1. Optimised RMS for a single aspheric surface at a given focal plane. 

#Compare with single spherical surface

#Initial guess for A4 and A6. Should be <<1 
A4 = 1e-6
A6 =0

#Set screen to be at 200 mm to compare with the single spherical case
focus_point = 200 


#Setting bound to look for solution. 
#Magnitude of A4 and A6 are definitely less than 1 
bound = ((-0.1,0.1), (-0.1,0.1))

#minimising
optimum = op.minimize(optimisation.aspheric_rms, [A4, A6], 
                      args=(100, 0.03, 1, 1.5,1/0.03, focus_point, 5), 
                      bounds=bound )

#Print out the optimised values 
optimum_param = optimum.x
min_rms = optimum.fun

print('The optimum A4, A6 is:', optimum_param)

print('The minimised RMS is:', min_rms, 'mm')

#%% 2. Optimised RMS for a single aspheric surface and planar surface.

#Compare with plano-convex lens 


#Initial guess for A4 and A6. Should be <<1 
A4 = 1e-6
A6 =0

#Set screen to be at 198.452 mm to compare with planoconvex
focus_point = 198.452700178 


#Setting bound to look for solution. Magnitude of A4 and A6 are definitely less than 1 
bound = ((-0.1,0.1), (-0.1,0.1))

#minimising
optimum = op.minimize(optimisation.plano_aspheric_rms, [A4, A6], 
                      args=(100, 5, 0.02, 1, 1.5168,1/0.02, focus_point ,10), 
                      bounds=bound )

optimum_param = optimum.x
min_rms = optimum.fun

print('The optimum (A4, A6) is:', optimum_param)

print('The minimised RMS is:', min_rms, 'mm')
        

#%% 3. Optimised RMS for a double aspheric surface lens.

#Compare with bi-convex lens 


#Initial guess for curvatures from biconvex optimisation

C1 = 0.01596323 #1st curvature
C2 = -0.00425303 #2nd curvature

#Initial guess for aspheric constants for 1st surface and 2nd surface.
A4  = 1e-9
A4_2 = -1e-9
A6  = A6_2 = 0


#Set screen to be at 198.452 mm to compare with bi-convex and plano-convex
focus_point = 198.452700178 


#Setting bound to look for solution. 
#Magnitude of A4 and A6 are definitely less than 1 
rmin = 5
bound = ( (0,1./rmin), (-1./rmin,0), (-0.1, 0.1),  
         (-0.1, 0.1), (-0.1, 0.1), (-0.1, 0.1) )

#minimising. Set maxfev high due to high number of parameters 
#Use Nelder-Mead method, which is not gradient based 

optimum = op.minimize(optimisation.double_aspheric_rms, [C1, C2, A4, 
                      A6, A4_2, A6_2], args=(100, 5, 1, 1.5168, focus_point, 10), 
                      bounds=bound, options = {'maxfev': 1e+3}, 
                      method = 'Nelder-Mead' )

#Print out the optimised values 

optimum_param = optimum.x
min_rms = optimum.fun

print('The optimum parameters are:', optimum_param)

print('The minimised RMS is:', min_rms, 'mm')

#10 times better than the initial bi-convex set up! Spot size 2e-4 mm ! 

