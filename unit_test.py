#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This script is used to perform unit test of the code in the raytracer.py file.

Run each cell one by one IN ORDER and make sure the output is as expected. 

Most cells: the user should observe the Console to see if 
                     the tests passed.

The first and last few cells: the user should examine the output and physics 
                            to see if it makes sense. The comments provide 
                            some guidances.

Created on Tue Oct 18 18:46:11 2022

@author: williamnguyen

"""
import raytracer as rt

import numpy as np  

import matplotlib.pyplot as plt


#%% Testing value error for class Ray () 

incorrect_position = [0, 0, 0 ,0]
direction = [1, 1, 1]
P1 = rt.Ray(incorrect_position, direction)

#check that a value error is raised in the console 

#%% Testing the str and repr method for class Ray () 

#Check the official and unofficial representation of Ray()
print('This unit test checks the official and unofficial representation of Ray()')

initial_position = [0,0,0]
initial_direction = [1,1,1]

P1 = rt.Ray(initial_position, initial_direction) 
official = repr(P1)
unofficial = str(P1)


print('\nOFFICIAL representation:', official)

print('\nUNOFFICIAL representation:', unofficial)

i = 0 

if repr(P1) == 'raytracer.Ray([0. 0. 0.],[1. 1. 1.])':
    print('\nOfficial representation is CORRECT')
    i +=1 
else: 
    print('\nOfficial representation is WRONG')
    

if str(P1) == 'Position vector is [0. 0. 0.] Direction vector is [1. 1. 1.]':
    print('\nUnofficial representation is CORRECT')
    i +=1 
else: 
    print('\nUnofficial representation is WRONG')
    
if i == 2:
    print('\nBoth representations are CORRECT')

#%% Testing p() and k() method of class Ray()
 
#Create an object of class Ray with given initial position and direction  
initial_position = [0,0,0]
initial_direction = [1,1,1]

P1 = rt.Ray(initial_position, initial_direction) 

print('Testing p() and k() of Ray class')

print('This unit test checks the position and direction returning method. \
      A Ray object is initialised with (0,0,0) and (1,1,1) as the position \
      and wave-vector. The methods return:')
      
print('\n INITIAL POSITION: ',P1.p())
print('\n WAVE VECTOR: ', P1.k())

print('\n RESULTS OF THE TESTS:')
if np.array_equal(P1.p(), initial_position) == True:
    print('\n Test 1 passed. p() method works!')
    
else: 
    print('\n Test 1 failed! p() method has bugs!')
    
if np.array_equal(P1.k(), initial_direction) == True:
    print('\n Test 2 passed. k() method works!')
    
else: 
    print('\n Test 2 failed! k() method has bugs!')
#%% Testing append() method of class Ray

#append new direction and point to the Ray object
new_point = [-1,-1,-1]
new_direction =  [1,2,3]
P1.append(new_point, new_direction) 

print('This unit test checks append() method of Ray class. A Ray object should\
      be appended with (-1,-1,-1) and (1,2,3) as the new position and wave-vector.')

print('The methods return:')    
      
print('\n NEW POSITION: ',P1.p())

print('\n NEW WAVE VECTOR: ', P1.k())

if np.array_equal(P1.k(), new_direction) == True and np.array_equal(
                                                    P1.p(), new_point) == True: 
    print('\n Test 3 passed. append() method works!')
else: 
    print('\n Test 3 failed! append() method has bugs!')
    
 

#%% Spherical lens interception testing (simple)

#Create 4 test Ray objects and a spherical lens object 


surface = rt.SphericalRefraction(z0 = 0, curvature = 1, n1 = 1, n2 = 1.5, 
                                 r_aperture = 0.5)

ray1 = rt.Ray([0,0,0], [0,0,1])
ray2 = rt.Ray([0,0,0], [1,1,0.15])
ray3 = rt.Ray([0,0,0], [1,1,0.1])
ray4 = rt.Ray([0,0,0],[-1,1,1]) 
 
intercept1 = surface.intercept(ray1) #should intercept
intercept2 = surface.intercept(ray2) #should intercept
intercept3 = surface.intercept(ray3) #should intercept
intercept4 = surface.intercept(ray4) #no intercept  (checked on geogebra)

print('This unit test checks intercept() method of SphericalRefraction class.')
print('A surface at z =0 with C = 1 and aperture radius 0.5 is used with 4 \
      test ray objects.')
print('The method returns:')
print("\n RAY 1 INTERCEPTION:", intercept1)
print("\n RAY 2 INTERCEPTION:", intercept2)
print("\n RAY 2 INTERCEPTION:", intercept3)
print("\n RAY 4 INTERCEPTION:", intercept4)
i = 0 
if len(intercept1) == 3:
    i+=1
    print('\n First intercept correct')
    
else: 
    print('\n First intercept wrong')   
    
if len(intercept2) == 3:
    i+=1 
    print('\n Second intercept correct')
    
    
else: 
    print('\n Second intercept wrong')   

if len(intercept3) == 3:
    i+=1 
    print('\n Third intercept correct')
   
else: 
    print('\n Third intercept wrong')   

if intercept4 == None:
    i+=1 
    print('\n Fourth intercept correct')
    
else: 
    print('\n Fourth intercept wrong')   

if i == 4:
    print('\n Test 4 passed! Simple intercept cases work!')


#%% Testing intercept with planar surface

planar_surface = rt.SphericalRefraction(100,0,1,1.5,50)
ray1 = rt.Ray([0,0,0], [0,0,1])
intercept = planar_surface.intercept(ray1)

print('This unit test checks the intercept method exception case.')
print('A ray starting at (0,0,0) moving +z should intercept a plane placed at \
      z0 = 100. Expected intercept is (0, 0, 100). The output is:')

print('\n INTERCEPT WITH PLANAR SURFACE IS ', intercept)
if np.array_equal(intercept, np.array([0, 0, 100])):
    print('\n Test 5 passed! Planar surface intercept works!')
else: 
    
    print('\n Test 5 failed! Planar surface has bugs!')

#%% Refraction Method Test
 
#Testing case: Incident ray hits at 45 degrees. 
#We expect refraction angle to  be 28 degrees.  

print('This unit test checks the refraction method via 3 cases.')
print('Case 1: Ray goes from air to glass at 45°. Should refract 28°. ')
print('Case 2: Ray goes from glass to air at 45°. Expect Total Internal \
      Reflection.')
print('Case 3: Ray goes from air to air 45°. No refraction expected.')
print('The ouputs were as follows:')

k = np.array([1,0,1])/np.sqrt(2)

n = np.array([0,0,-1]) 

out_1 = rt.refract(k, n, 1, 1.5) #air to glass
out_2 = rt.refract(k, n, 1.5, 1) #glass to air
out_3 = rt.refract(k, n, 1, 1) #air to air 

print('\n CASE 1 REFRACTED RAY DIRECTION : ', out_1)
print('\n CASE 2 REFRACTED RAY DIRECTION : ', out_2)
print('\n CASE 3 REFRACTED RAY DIRECTION : ', out_3)


refracted_angle_1 = np.arccos(np.abs(np.dot(out_1, n))) /np.pi * 180
refracted_angle_3 = np.arccos(np.abs(np.dot(out_3, n))) /np.pi * 180


print('\n CASE 1 ANGLE OF REFRACTION (°): ',refracted_angle_1)

if out_2 == None:
    print('\n CASE 2: Total internal reflection')
else:
    print('\n CASE 2: NO total internal reflection. Bugs!')
    
print('\n CASE 3 ANGLE OF REFRACTION (°): ',refracted_angle_3)

if (np.round(refracted_angle_1,0) == 28 and refracted_angle_3 == 45.0 
    and out_2 == None):
    
    print('\n Test 6 passed! refract() method works!')
    
else: 
    print('\n Test 6 failed! refract() method has bugs!')

#%% Testing Bundle Class. 


bundle1 = rt.Bundle(5,10).make_bundle() #should be a ray object

#returns position of the first ray in bundle. Should be (-2.5, 0, 0)
print('EXPECTED INITIAL POSITION OF FIRST RAY IN THE BUNDLE IS',(-2.5, 0, 0))

print('INITIAL POSITION OF FIRST RAY IN THE BUNDLE IS',bundle1[0].p()) 

#%% Testing the WARNING from the Bundle Class. 
 
bundle3 = rt.Bundle(5,123).make_bundle()


#warning can be seen on IDE. Works! 

#%% Tracing a few rays through surface to test model (Task 9)

"""
Getting started: 
    - Input any value of non-zero value of C to see how the model behaves 
    - Observe the graph produced: where do the rays intercept? 
                                  do they rays converge or diverge? 
    - To test, use C = 0.03 and C = -0.03 
    
For C = 0.03, we expect: 
    
    - Rays converge at 200 mm
    - Rays intercept the convex curve correctly

For C = -0.03, we expect

    - Rays diverge after hitting the concave surface

"""
C = 0.03 #input 0.03 or -0.03 to test. Can play with other values
R = 1/np.abs(C)

surface = rt.SphericalRefraction(z0=100, curvature = C, 
                                 n1= 1, n2=1.5, r_aperture = R)

screen = rt.OutputPlane(z0 = 200)

ray1 = rt.Ray([-5,0,0], [0,0,1])
ray2 = rt.Ray([5,0,0], [0,0,1])
ray3 = rt.Ray([10,0,0], [0,0,1])
ray4 = rt.Ray([-10,0,0], [0,0,1])

ray = (ray1, ray2, ray3, ray4)
elements = (surface, screen)


for i in range(len(ray)):
    for elem in elements:
        elem.propagate_ray(ray[i])
        
        
x1, z1 = [], []

for j in range(len(ray1.verticies())):
    x1.append(ray1.verticies()[j][0])
    z1.append(ray1.verticies()[j][2])
    
x2, z2 = [], [] 
for j in range(len(ray2.verticies())):
    x2.append(ray2.verticies()[j][0])
    z2.append(ray2.verticies()[j][2])

x3, z3 = [], [] 
for j in range(len(ray3.verticies())):
    x3.append(ray3.verticies()[j][0])
    z3.append(ray3.verticies()[j][2])
    
x4, z4 = [], [] 
for j in range(len(ray4.verticies())):
    x4.append(ray4.verticies()[j][0])
    z4.append(ray4.verticies()[j][2])
    

plt.plot(z1,x1, marker='x')
plt.plot(z2,x2, marker='x')
plt.plot(z4,x4, marker='x')
plt.plot(z3,x3, marker='x')


plt.xlabel('z / mm', fontsize=15)
plt.ylabel('x / mm ',fontsize=15)

plt.axvline(x = 200, color = 'black', label = 'Screen')
zplot = np.linspace(-R, R,10000)
if surface._curvature > 0:
    plt.plot(100 + R -np.sqrt(R**2 - zplot**2),zplot,label='Surface')

if surface._curvature < 0:
    plt.plot(100 - R +np.sqrt(R**2 - zplot**2),zplot,label='Surface')


plt.legend()

 

#%% Testing intercept between 2D rays to find paraxial focal plane

"""
This tests a simple code used to calculate intercept between 2D rays. This
is used to find the paraxial focal plane. 

The two rays were used. The expected z intercept is 1.0 (solved on paper). 

"""


def slope(point1, point2):
     return (point2[1]- point1[1])/(point2[0] - point1[0])

def yintercept(point1, point2):
     return point2[1] - slope(point1, point2) * point2[0]

def intercept(ray1, ray2): #ray 1 and ray 2 are Ray objects. 
    #Returns z coordinate of intercept
     #setting up points

     point1_ray1 = np.array( [ray1.p()[2] , ray1.p()[0] ] )
     point2_ray1 = ray1.p() + ray1.k()
     point2_ray1  = np.array( [point2_ray1 [2] , point2_ray1 [0] ])

     point1_ray2 = np.array( [ray2.p()[2] , ray2.p()[0] ] )
     point2_ray2 = ray2.p() + ray2.k()
     point2_ray2  = np.array( [point2_ray2 [2] , point2_ray2 [0]] )
     
     slope1 = slope(point1_ray1, point2_ray1 )
     slope2 = slope(point1_ray2, point2_ray2 )
     
     if slope1 == slope2: 
         print('Parallel rays do not intersect')
         return None
     
     else:   
         xintercept1 = yintercept(point1_ray1, point2_ray1 )
         xintercept2 = yintercept(point1_ray2, point2_ray2 )
         
         # z coordinate of intercept
         z_intercept = (xintercept2 - xintercept1)/(slope1 - slope2)
         
         return z_intercept
 
 
ray1 = rt.Ray(np.array([0.1, 0, 0]), np.array([0,0,1]))
ray2 = rt.Ray(np.array([0, 0, 0]), np.array([0.1,0,1]))
z_intercept = intercept(ray1, ray2)
print('z intercept coordinate', z_intercept)

if z_intercept == 1.:
    print('Test passed! Method to find intercept between rays works!')
else:
    print('Test failed! Has bugs!')
    
#%% Testing automatically finding paraxial focal plane.
"""
This tests the implementation of the functions used in the cell above to find
the paraxial focal plane automatically.

The test set up is as follows: 
    1. Convex surface placed at z = 100 
    2. Concave surface placed at z = 105  
    3. Screen placed at 170 

Two paraxial rays are traced through the optical system and their intercept 
is found as the focal point. 

The result is plotted. Observe and check the intercept is correct. 

"""

#Setting up optical elements

n1 = 1
n2 = 1.5168
R1 = 1/0.03
R2 = 1/0.03

#Initialise the optical elements 
surface_1 = rt.SphericalRefraction(100, 1/R1, n1, n2, R1)
surface_2 = rt.SphericalRefraction(105, -1/R2, n2, n1, R2)
screen = rt.OutputPlane(170)
elements = (surface_1, surface_2, screen)
 
#Set up paraxial rays
ray1 = rt.Ray(np.array([1,0,0]), np.array([0,0,1]))
ray2 = rt.Ray(np.array([-1,0,0]), np.array([0,0,1]))

#Propagating through the optical elements 
for elem in elements:
     elem.propagate_ray(ray1)
     elem.propagate_ray(ray2)
 
#find intercept between two paraxial rays
z_intercept = rt.intercept(ray1, ray2)

x1, z1 =[], []

for j in range(len(ray1.verticies())):
    x1.append(ray1.verticies()[j][0])
    z1.append(ray1.verticies()[j][2])
    
x2, z2 = [], []
for j in range(len(ray2.verticies())):
    x2.append(ray2.verticies()[j][0])
    z2.append(ray2.verticies()[j][2])

#Plotting to see that everything makes sense 

plt.plot(z1,x1,marker='x')
plt.plot(z2,x2,marker='x')


zplot = np.linspace(-33.3, 33.3, 1000)
plt.plot(100 + R1 -np.sqrt(R1**2 - zplot**2),zplot, label ='Surface 1')
plt.plot(105 - R2 +np.sqrt(R2**2 - zplot**2),zplot, label ='Surface 2')
plt.axvline(z_intercept, label='Paraxial focal plane') 
#make sure this is where the paraxial ray intersects

plt.legend()
plt.xlabel('z / mm')
plt.ylabel('x / mm')    

#You can see that the intercept function to find the paraxial focus works 


    