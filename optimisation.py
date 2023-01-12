#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains the optimisation methods to minimise the RMS of different
optical systems. 

Created on Sat Dec  3 23:45:24 2022

@author: williamnguyen
"""


import numpy as np
import raytracer as rt 
    
#Lens optimisation function

def rms(x, y):
    """
    FUNCTION THAT CALCULATES RMS SPOT SIZE GIVEN THE LIST OF X AND Y POSTIONS

    Parameters
    ----------
    x : x position of the rays at the chosen plane 
    (list)
    y : y position of the rays at the chosen plane 
    (list)   .
    
    
    Returns
    -------
    (rms spot size, uncertainty)

    """
    x = np.array(x)
    y = np.array(y)
    distance = np.sqrt(x**2 + y**2)
    rms_val = np.sqrt(np.mean(distance**2) )
    
    
    return rms_val



def system_rms(c, z0, z1, z2, D):
    
    """ FUNCTION THAT CALCULATES RMS AT A GIVEN POINT 
    IN AN OPTICAL SYSTEM (glass and air)
    
    Parameters: 
    ----------
    
    C1: curvature of surface 1 
    C2: curvature of surface 2
    z0: position of surface 1 (on the z axis)
    z1: position of surface 2(on the z axis)
    z2: position of the screen (on the z axis)
    D: diameter of the beam bundle 
    
    Returns:
    ----------
    
    RMS (geometrical spot size) at the screen  
    
            """
    C1 = c[0]
    C2 = c[1]        
    bundle = rt.Bundle(D,10).make_bundle()
    
    if C1 == 0: 
        surface1 = rt.SphericalRefraction(z0, C1, 1, 1.5168, 60)
    else: 
        surface1 = rt.SphericalRefraction(z0, C1, 1, 1.5168, np.abs(1/C1))
    
    if C2 == 0: 
        surface2 = rt.SphericalRefraction(z1, C2, 1.5168, 1, 60)
    else: 
        surface2 = rt.SphericalRefraction(z1, C2, 1.5168, 1, np.abs(1/C2))
        
    screen = rt.OutputPlane(z2)
    
    elements = (surface1, surface2, screen)

    for i in range(len(bundle)):
        for elem in elements:
            elem.propagate_ray(bundle[i])
            
    x = []
    y = []

    for k in range(len(bundle)):
        if bundle[k]._terminate == False: 
            x.append(bundle[k].p()[0])
            y.append(bundle[k].p()[1])
    
    return rms(x,y)
    

def aspheric_rms(A, z0, C, n1, n2, r_aperture, screen_position, D = 10):
    """ 
    FUNCTION THAT CALCULATES RMS SPOT SIZE FOR 
    A GIVEN ASPHERIC SURFACE AND SCREEN DISTANCE
    
    Input
    -----
    screen_position: z position of the screen. 
    
    diameter: diameter of the beam. Default is 10 mm
    
    Output
    ------
    
    rms spot size associated with the apsheric surface at a given distance. 
    
    """    
    A_4 = A[0]
    A_6 = A[1]
    
    bundle = rt.Bundle(D,10).make_bundle()
 
    surface = rt.AsphericRefraction(z0, C, A_4, A_6, n1, n2, r_aperture)       
    screen = rt.OutputPlane(screen_position)
    
    elements = (surface, screen)

    for i in range(len(bundle)):
        for elem in elements:
            elem.propagate_ray(bundle[i])
            
    x = []
    y = []

    for k in range(len(bundle)):
        if bundle[k]._terminate == False: 
            x.append(bundle[k].p()[0])
            y.append(bundle[k].p()[1])
    
    return rms(x,y)


def plano_aspheric_rms(A, z0, sep, C, n1, n2, r_aperture, 
                       screen_position, D = 10):
    """ 
    FUNCTION THAT CALCULATES RMS SPOT SIZE FOR A GIVEN 
    PLANO-ASPHERIC LENS AT A GIVEN SCREEN POSITION 
    
    Input
    -----
    screen_position: z position of the screen. 
    
    sep = distance between planar surface and aspheric surface 
    
    diameter: diameter of the beam. Default is 10 mm
    
    Output
    ------
    
    rms spot size associated with the apsheric surface at a given distance. 
    
    """    
    A_4 = A[0]
    A_6 = A[1]
    
    bundle = rt.Bundle(D,10).make_bundle()
    
    
    surface = rt.AsphericRefraction(z0, C, A_4, A_6, n1, n2, r_aperture) 
    plane = rt.SphericalRefraction(z0 + sep, 0, n2, n1, r_aperture)      
    screen = rt.OutputPlane(screen_position)
    
    elements = (surface, plane, screen)

    for i in range(len(bundle)):
        for elem in elements:
            elem.propagate_ray(bundle[i])
            
    x = []
    y = []

    for k in range(len(bundle)):
        if bundle[k]._terminate == False: 
            x.append(bundle[k].p()[0])
            y.append(bundle[k].p()[1])
    
    return rms(x,y)


def aspheric_sphere_rms(A, z0, sep, C1, C2, n1, n2, r_aperture, 
                        screen_position, D = 10):
    
    A_4 = A[0]
    A_6 = A[1]
    
    bundle = rt.Bundle(D,10).make_bundle()
    
    
    surface = rt.AsphericRefraction(z0, C1, A_4, A_6, n1, n2, r_aperture) 
    surface_sphere = rt.SphericalRefraction(z0 + sep, C2, n2, n1, np.abs(1/C2))      
    screen = rt.OutputPlane(screen_position)
    
    elements = (surface, surface_sphere, screen)

    for i in range(len(bundle)):
        for elem in elements:
            elem.propagate_ray(bundle[i])
            
    x = []
    y = []

    for k in range(len(bundle)):
        if bundle[k]._terminate == False: 
            x.append(bundle[k].p()[0])
            y.append(bundle[k].p()[1])
    
    return rms(x,y)

    

def double_aspheric_rms(param, z0, sep, n1, n2, screen_position, D = 10):
    """ 
    FUNCTION THAT CALCULATES RMS SPOT SIZE FOR A GIVEN DOUBLE ASPHERIC SURFACE 
    LENS AT A GIVEN SCREEN POSITION 
    
    Input
    -----
    
    param: array in the form of [C1, C2, A_4, A_6, A_4_2, A_6_2].
           Parameters are to be optimised to minimise rms spot size. 
    
    
    screen_position: z position of the screen. 
    
    sep = distance between planar surface and aspheric surface 
    
    diameter: diameter of the beam. Default is 10 mm
    
    Output
    ------
    
    rms spot size associated with the apsheric surface at a given distance. 
    
    """    
    #curvatures of two surfaces
    C1 = param[0] 
    C2 = param[1]
    
    #A4 and A6 of two surfaces 
    A_4 = param[2]
    A_6 = param[3]
    A_4_2 = param[4] 
    A_6_2 = param[5] 
    
    
    
    if C1 == 0: 
        surface = rt.AsphericRefraction(z0, C1, A_4, A_6, n1, n2, 60) 
    
    else: 
        surface = rt.AsphericRefraction(z0, C1, A_4, A_6, n1, n2, np.abs(1/C1)) 

    
    if C2 == 0: 
        surface2 = rt.AsphericRefraction(z0 + sep, C2, A_4_2, A_6_2, n2, n1, 60)  
        
    else: 
        surface2 = rt.AsphericRefraction(z0 + sep, C2, A_4_2, A_6_2, 
                                      n2, n1, np.abs(1/C2))
        
    screen = rt.OutputPlane(screen_position)
    
    elements = (surface, surface2, screen)
    
    bundle = rt.Bundle(D,10).make_bundle()
    
    for i in range(len(bundle)):
        for elem in elements:
            elem.propagate_ray(bundle[i])
            
    x = []
    y = []

    for k in range(len(bundle)):
        if bundle[k]._terminate == False: 
            x.append(bundle[k].p()[0])
            y.append(bundle[k].p()[1])
    
    
    return rms(x,y)
