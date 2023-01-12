#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file contains the code to investigate ray tracing: 
    
- Tracing ray bundle through a single spherical surface
- Comparing the two different orientations of a plano-convex lens

Remember running the cell in different orders will affect the outputs. 

Created on Mon Oct 31 10:56:42 2022

@author: williamnguyen.
"""

import raytracer as rt
import optimisation 
import numpy as np 
import matplotlib.pyplot as plt

#%% Tracing the bundle . Task 12

#A single spherical surface, beyond the paraxial limit

#Set up 
surface = rt.SphericalRefraction(z0 = 100, curvature = 0.03, n1 = 1, 
                                 n2 = 1.5, r_aperture = 1/0.03)

screen = rt.OutputPlane(z0 = 200) #placed at paraxial focus 

bundle1 = rt.Bundle(5, 10).make_bundle()  #5 mm beam
elements = (surface, screen) 

#Propagate the beam through the optical elements 

for i in range(len(bundle1)):
    for elem in elements:
        elem.propagate_ray(bundle1[i])


 
for j in range(len(bundle1)):
    x=[ bundle1[j].verticies()[0][0], bundle1[j].verticies()[1][0], 
       bundle1[j].verticies()[2][0] ]
    
    z = [ bundle1[j].verticies()[0][2], bundle1[j].verticies()[1][2], 
         bundle1[j].verticies()[2][2] ]
    
    plt.plot(z,x,color='blue')
    plt.xlabel('z / mm',fontsize=15)
    plt.ylabel('x / mm',fontsize=15)

#Plotting optical elements
plt.axvline(x=200, label='Screen',color='black')
zplot = np.linspace(-33.3,33.3,10000)
plt.plot(133.33-np.sqrt(33.3**2 - zplot**2),zplot,label='Surface')
plt.ylim(-10,10)
plt.legend()
plt.savefig('Task12.png',dpi=500, bbox_inches='tight')

#%% Spot diagram at the focal plane ( z = 200 ). Task 13

x =[]
y = []

for k in range(len(bundle1)):
    x.append(bundle1[k].p()[0])
    y.append(bundle1[k].p()[1])

plt.scatter(x,y,color='blue', marker='o')
plt.xlabel('x / mm',fontsize=15)
plt.ylabel('y / mm',fontsize=15)

rms_1 = optimisation.rms(x,y)
plt.savefig('spotdiagram_focalpoint.png',dpi=500, bbox_inches='tight')
print('The rms value is', rms_1, 'mm')


#%% Spot diagram at the starting  point ( z = 0 ). Task 13

x =[]
y = []

for k in range(len(bundle1)):
    x.append(bundle1[k].verticies()[0][0])
    y.append(bundle1[k].verticies()[0][1])

plt.scatter(x,y,color='red', marker='o')
plt.xlabel('x / mm',fontsize=15)
plt.ylabel('y / mm',fontsize=15)
#plt.title('Spot diagram at z = 0')
plt.savefig('spotdiagram_z0.png',dpi=500, bbox_inches='tight')


#%% Modelling a plano-convex singlet lens. Task 15. 


#set beam diameter. Investigate 1 mm - 10 mm, increment 1 mm

d = 10

#set up. First orientation: planar surface followed by curved surface

plane = rt.SphericalRefraction(z0 = 100, curvature = 0, n1= 1, 
                               n2 = 1.5168, r_aperture =1/0.02)

convex = rt.SphericalRefraction(z0 = 105, curvature= -0.02, n1= 1.5168,
                                n2 = 1, r_aperture =1/0.02)

screen = rt.OutputPlane(201.7487808269458) #at paraxial focus

elements = (plane, convex, screen)

bundle2 = rt.Bundle(d,10).make_bundle() 

#Propagate the rays through the optical elements 

for i in range(len(bundle2)):
    for elem in elements:
        elem.propagate_ray(bundle2[i])


for j in range(len(bundle2)):
    
    x=[ bundle2[j].verticies()[0][0], bundle2[j].verticies()[1][0], 
       bundle2[j].verticies()[2][0], bundle2[j].verticies()[3][0] ]
    
    z = [ bundle2[j].verticies()[0][2], bundle2[j].verticies()[1][2], 
         bundle2[j].verticies()[2][2], bundle2[j].verticies()[3][2] ]
    
    plt.plot(z,x,color='blue')
    plt.xlabel('z / mm',fontsize=15)
    plt.ylabel('x / mm',fontsize=15)

#plotting the surfaces
zplot = np.linspace(-50, 50, 1000)
plt.plot(105 - 50 + np.sqrt(50**2 - zplot**2), zplot, label = 'curved surface')
plt.axvline(x=100, ymin = -40, ymax=40, label = 'planar surface', color = 'red')
plt.axvline(x=201.74878082694576, ymin = -40, ymax=40, 
            label = 'screen at z $\sim 202$', color = 'black')
plt.ylim(-20,20)
plt.legend()

#Finding the paraxial focus 
ray1 = rt.Ray(np.array([0.1,0,0]), np.array([0,0,1]) ) 
ray2 = rt.Ray(np.array([-0.1,0,0]), np.array([0,0,1]) ) 

for elem in elements:
    elem.propagate_ray(ray1)
    elem.propagate_ray(ray2)

paraxial_focus = rt.intercept(ray1, ray2)

print('the paraxial focus is at', paraxial_focus)
#%% 
# Spot diagram at z = 201.7487808269458 mm. First configuration

x =[]
y = []

for k in range(len(bundle2)):
    x.append(bundle2[k].verticies()[3][0])
    y.append(bundle2[k].verticies()[3][1])

plt.scatter(x,y,color='blue', marker='o')
plt.xlabel('x / mm',fontsize=15)
plt.ylabel('y / mm',fontsize=15)
plt.title('Spot diagram at the focal plane z = 201.7 mm. First configuration')

rms = optimisation.rms(x, y)

plt.savefig('first_configuration.png',dpi=500)

print('The rms value is', rms, 'mm')

#%% Curved - Plane configuration (Swapped position of plane and convex )

convex = rt.SphericalRefraction(z0=100, curvature=0.02, n1= 1, 
                                n2=1.5168, r_aperture =1/0.02)

plane = rt.SphericalRefraction(z0=105, curvature=0, n1= 1.5168, 
                               n2=1, r_aperture =1/0.02)

 
bundle2 = rt.Bundle(d,10).make_bundle() 
screen = rt.OutputPlane(198.45270017751582) #placed at paraxial focus
elements = (convex, plane, screen)
#Propagate the rays through the optical elements 

for i in range(len(bundle2)):
    for elem in elements:
        elem.propagate_ray(bundle2[i])


for j in range(len(bundle2)):
    
    x=[ bundle2[j].verticies()[0][0], bundle2[j].verticies()[1][0], 
       bundle2[j].verticies()[2][0], bundle2[j].verticies()[3][0] ]
    
    z = [ bundle2[j].verticies()[0][2], bundle2[j].verticies()[1][2], 
         bundle2[j].verticies()[2][2], bundle2[j].verticies()[3][2] ]
    
    plt.plot(z,x,color='blue')
    plt.xlabel('z / mm',fontsize=15)
    plt.ylabel('x / mm',fontsize=15)
    
#%% Spot diagram at the focal plane ( z = 198.45270017751582 mm ). 
x =[]
y = []

for k in range(len(bundle2)):
    x.append(bundle2[k].verticies()[3][0])
    y.append(bundle2[k].verticies()[3][1])

plt.scatter(x,y,color='blue', marker='o')
plt.xlabel('x / mm',fontsize=15)
plt.ylabel('y / mm',fontsize=15)
plt.title('Spot diagram at the focal plane. Second configuration')
x = np.array(x)
y = np.array(y)
rms_2 = np.sqrt(np.mean(x**2 + y**2))
plt.savefig('second_configuration.png',dpi=500)
print('The rms value is', rms_2, 'mm')

#%% Comparing the performances of the two configruations (consider RMS)

#Remember to use your own filepath! 

diameter, rms1, rms2 = np.loadtxt('/Users/williamnguyen/Desktop/Y2/optical ray tracer/lens_performance.txt',
                                  unpack=True, skiprows=1)

plt.errorbar(diameter, rms1*1000, yerr = rms1*100,
             label='Plane-surface configuration', fmt='.k', color='blue',
             capsize = 3)

plt.errorbar(diameter, rms2*1000, yerr = rms2*100,
             label='Surface-plane configruation',fmt='.k', color='red',
             capsize = 3)

plt.ylabel('RMS spot size / $\mu m$', fontsize=13)
plt.xlabel('Beam diameter / mm', fontsize=13)
plt.grid()
plt.legend()
plt.savefig('performance.png',dpi=500, bbox_inches='tight')

#Uncertainty of 10% for RMS value. 100 rays is about 10% deviation from when
#a million rays are used. 
