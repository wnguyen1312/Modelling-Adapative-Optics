#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This module contains optimisation and testing for the aspheric surface: 
        - Reduce the apsherical class to a spherical class by setting 
        A_4 and A_6 to 0.
        - Check that the output is the same as before.
Created on Thu Nov 17 20:04:28 2022

@author: williamnguyen
"""

import raytracer as rt
import numpy as np 
import matplotlib.pyplot as plt
import optimisation 

#We expect focal point to occur at z = 200 mm 

A4 = 0 
A6 = 0
surface = rt.AsphericRefraction(100, 0.03, A4, A6, 1, 1.5, 1/0.03)

#surface2 = rt.AsphericRefraction(105, -0.03,A4,A6,1.5,1,1/0.03)


screen = rt.OutputPlane(z0 = 200)

elements = (surface, screen)

ray1 = rt.Ray(np.array([-1,0,0]),np.array([0,0,1]))
ray2 = rt.Ray(np.array([1,0,0]),np.array([0,0,1]))
ray3 = rt.Ray(np.array([-2,0,0]),np.array([0,0,1]))
ray4 = rt.Ray(np.array([2,0,0]),np.array([0,0,1]))

ray = (ray1, ray2, ray3, ray4)


for i in range(len(ray)):
    for elem in elements:
        elem.propagate_ray(ray[i])
        
x1, y1, z1 =[], [], []

for j in range(len(ray1.verticies())):
    x1.append(ray1.verticies()[j][0])
    y1.append(ray1.verticies()[j][1])
    z1.append(ray1.verticies()[j][2])
    
x2, y2, z2 = [], [], []
for j in range(len(ray2.verticies())):
    x2.append(ray2.verticies()[j][0])
    y2.append(ray2.verticies()[j][1])
    z2.append(ray2.verticies()[j][2])

x3, y3, z3 = [], [], []
for j in range(len(ray3.verticies())):
    x3.append(ray3.verticies()[j][0])
    y3.append(ray3.verticies()[j][1])
    z3.append(ray3.verticies()[j][2])
    
x4, y4, z4 = [], [], []
for j in range(len(ray4.verticies())):
    x4.append(ray4.verticies()[j][0])
    y4.append(ray4.verticies()[j][1])
    z4.append(ray4.verticies()[j][2])
    
fig = plt.figure(figsize=(6.4, 4.8))
plt.plot(z1,x1, marker='x')
plt.plot(z2,x2, marker='x')
plt.plot(z4,x4, marker='x')
plt.plot(z3,x3, marker='x')


plt.xlabel('z / mm', fontsize=15)
plt.ylabel('x / mm ',fontsize=15)

zplot = np.linspace(-33.3,33.3,1000)

compare = 100 + 0.03*zplot**2/(1+np.sqrt(1-0.03*0.03*zplot*zplot)) 
+ A4*zplot**4 + A6*zplot**6

plt.plot(compare,zplot)


#%% Tracing 5 mm beam through the apsheric surface. Compare with spherical case.

bundle1 = rt.Bundle(5, 10).make_bundle() 

#Propagate the rays through the optical elements 

for i in range(len(bundle1)):
    for elem in elements:
        elem.propagate_ray(bundle1[i])


 
for j in range(len(bundle1)):
    
    x= [ bundle1[j].verticies()[0][0], bundle1[j].verticies()[1][0], 
       bundle1[j].verticies()[2][0] ]
    
    z = [ bundle1[j].verticies()[0][2], bundle1[j].verticies()[1][2], 
         bundle1[j].verticies()[2][2] ]
    
    plt.plot(z,x,color='blue')
    plt.xlabel('z / mm',fontsize=15)
    plt.ylabel('x / mm',fontsize=15)
    
plt.plot(compare,zplot)
plt.ylim(-4,4)
plt.axvline(x=200, label='Screen',color='black')

#%% Spot diagram at the focal plane ( 200 mm)

x =[]
y = []

for k in range(len(bundle1)):
    x.append(bundle1[k].p()[0])
    y.append(bundle1[k].p()[1])

plt.scatter(x,y,color='blue', marker='o')
plt.xlabel('x / mm',fontsize=15)
plt.ylabel('y / mm',fontsize=15)
plt.title('Spot diagram at the focal plane z = 200')

rms_1 = optimisation.rms(x,y)
print('The rms value is', rms_1, 'mm')

#We expect 0.0015685465095563694 mm (as per the spherical surface case)


