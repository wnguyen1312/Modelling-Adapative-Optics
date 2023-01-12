#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This file contains the instruction on how to get started with raytracing and
visualising the results: 
    
1. How to get set up an optical system
2. How to set up rays or beams
3. How to propagate rays or beams 
4. Automatically finding the paraxial focus 
4. Visualising the trajectory 
5. Spot diagram and rms value at the focal plane

See raytracer.py for more info on each class. 

Created on Wed Nov 23 18:54:47 2022

@author: williamnguyen
"""

#import important modules

import raytracer as rt
import numpy as np 
import matplotlib.pyplot as plt
import optimisation

#%% 1. How to set up an optical system 

"""
There are 3 types of optical elements to initialise: 
    - Spherical surface: from the SphericalRefraction class
    - Aspheric surface:  from the AsphericRefraction class
    - Output plane: from the OutputPlane class
"""

#Example of each optical element: 
    
#convex surface at z = 100 mm 
convex_surface = rt.SphericalRefraction(z0 = 100, curvature = 0.03, n1 = 1, 
                                           n2 = 1.5, r_aperture = 1/0.03)

#concave surface at z = 105 mm 
concave_surface = rt.SphericalRefraction(z0 = 105, curvature = -0.03, n1 = 1.5, 
                                           n2 = 1, r_aperture = 1/0.03)

#aspheric surface at z = 100 mm 
aspheric_surface = rt.AsphericRefraction(z0 = 105, curvature = -0.03, A_4 = 1e-6, 
                                          A_6 = -1e-8, n1 = 1, 
                                          n2 = 1.5, r_aperture = 1/0.03)
#screen at z = 137.43325742809589 mm 
screen = rt.OutputPlane(z0 = 137.43325742809589)

#Put everything into a list to make an optical system
#Here we will a convex and concave surface to make a bi-convex lens

elements = (convex_surface, concave_surface, screen)
#notice we didn't include the aspheric_surface for this example. 

#%% 2. Setting up your rays or beams 

#Paraxial Ray. 
#Initial position at (±1, 0, 0) moving in the + z direction. 

ray_1 = rt.Ray(r_position = np.array([1, 0, 0]), 
               r_direction = np.array([0, 0, 1]))
ray_2 = rt.Ray(r_position = np.array([-1, 0, 0]), 
               r_direction = np.array([0, 0, 1]))

#Beam. 5 mm diameter beam consisting of 10^2 individual rays. 

beam = rt.Bundle(diameter = 5, n = 10).make_bundle()

#%% 3. Propagating your rays or your beam

#Propagate the ray/beam through each element 1 by 1 

#Ray
for elem in elements:
    elem.propagate_ray(ray_1)
    
for elem in elements:
    elem.propagate_ray(ray_2)

#Beam  

for i in range(len(beam)):
    for elem in elements:
        elem.propagate_ray(beam[i])

#%% 4. Finding the paraxial focus 

#The paraxial focus is the intercept between two paraxial rays

paraxial_focus = rt.intercept(ray_1 , ray_2)

print('The focus is at', paraxial_focus)
#%% 5. Visualising the trajectories. We will do it for a 5 mm beam 

for j in range(len(beam)):
    x=[ beam[j].verticies()[0][0], beam[j].verticies()[1][0], 
       beam[j].verticies()[2][0], beam[j].verticies()[3][0] ]
    
    z = [ beam[j].verticies()[0][2], beam[j].verticies()[1][2], 
         beam[j].verticies()[2][2], beam[j].verticies()[3][2] ]
    
    plt.plot(z,x,color='blue')
    plt.xlabel('z / mm',fontsize=15)
    plt.ylabel('x / mm',fontsize=15)

#Plotting optical elements for visualisation 
plt.axvline(x = paraxial_focus, label='Screen',color='black')
zplot = np.linspace(-20,20,10000)
plt.plot(133.33-np.sqrt(33.3**2 - zplot**2),zplot,label='First surface')
plt.plot(105-33.333+np.sqrt(33.3**2 - zplot**2),zplot,label='Second surface')
plt.legend()

#%% 5. Getting the spot size and spot diagram at the focal plane

x =[]
y = []

for i in range(len(beam)):
    x.append(beam[i].p()[0])
    y.append(beam[i].p()[1])

plt.scatter(x,y,color='blue', marker='o')
plt.xlabel('x / mm',fontsize=15)
plt.ylabel('y / mm',fontsize=15)

rms_1 = optimisation.rms(x,y)

print('RMS spot size is', rms_1 , 'mm')