# -*- coding: utf-8 -*-
"""

This module contains all the essential classes and functions for ray tracing.

Created on Mon Oct 17 09:31:50 2022

@author: William Nguyen 
"""

import numpy as np
import scipy.optimize as op 
import warnings 

class Ray:
    """
    CLASS USED TO REPRESENT A LIGHT RAY
    
    Attributes: 
    ----------
    
    r_position (list len 3): initial position vector of the ray
    
    r_direction (list len 3): direction vector of the ray
    
    terminate (boolean): if True then ray is terminated (no more propagation)
    
    Methods: 
    ----------  
    
    p: returns the ray's current position vector
    k: returns the ray's current direction vector
    append: store the new point and direction to an array
    verticies: returns all the points on the ray 
    
    """
    
    def __init__(self, r_position, r_direction):
        
        if len(r_position) != 3 or len(r_position) != 3 :
            raise ValueError("Position & direction list needs a length of 3.")
            
        self.__position = np.array(r_position, dtype = float)
        self.__dir = np.array(r_direction, dtype = float)
        self.__points = [self.__position]
        self.__direction = [self.__dir]
        self._terminate = False 
        
    def __str__(self): 
        
        return('Position vector is ' + str(self.__points[0]) + 
               ' Direction vector is ' + str(self.__direction[0]) )
    
    def __repr__(self):
        
        return('raytracer.Ray(' + str(self.__points[0]) + ','
               + str(self.__direction[0]) + ')')
        
    def p(self):
        
        return self.__points[-1] 
    
    def k(self):
        
        return self.__direction[-1]
    
    def append(self, p, k):
        
        self.__points.append(p)
        self.__direction.append(k)
        
    def verticies(self):
        
        return self.__points



#Optical element class used to represent surfaces, lenses, etc

#Base class 
class OpticalElement:
    
    def propagate_ray(self, ray):
        "propagate  ray through an optical element"
        raise NotImplementedError()



#Refraction function.
def refract(k_in, n, n1, n2):
    """ 
    FUNCTION CALCULATES THE REFRACTED RAY DIRECTION BY SNELL'S LAW. 
    RAY TRAVELS FROM MEDIUM 1 --> MEDIUM 2. 
    
    Parameters: 
    ----------
    k_in = normalised incident wave-vector 
    n = normalised surface normal pointing OUT of the surface 
    n1 = refractive index of medium 1
    n2 = refractive index of medium 2 
    
    theta_1 = incident angle found via dot product 
    k_out= normalised refracted wave-vector

    Returns:
    ----------
    
    For a given incident angle, returns: 
        k_out
    If no valid solution, returns:
        None
    """
    
    k_in, n = np.array(k_in, dtype = float), np.array(n, dtype = float)
    
    #Find the incident angle to check for total internal reflection 
    if len(k_in) != 3: 
        raise TypeError("Wave-vector input should have a length of 3")
    
    theta_1 = np.arccos(np.abs(np.dot(n,k_in))) 
    
    if np.sin(theta_1) > n2/n1: 
        print('Total internal reflection occurs')
        return None 
    
    else: #Refraction
        u = n1 / n2 #
        term1 = u*(np.cross(n, np.cross(-n,k_in)))
        term2 = -n*np.sqrt(1-u*u *np.dot(np.cross(n,k_in), np.cross(n,k_in)) )
        k_out = term1 + term2 #refracted wave-vector
        return k_out 

#Function to check Line & Plane Collision to be used in Spherical Refraction

def line_plane_intercept(plane_norm, plane_point, ray_direction, 
                         ray_point):
    """
    FUNCTION USED TO CHECK FOR INTERCEPT BETWEEN A LINE AND A PLANE. TO BE USED
    FOR SPECIAL CASE OF A PLANAR SURFACE IN SPHERICAL REFRACTION CLASS. 
    
    Input: 
    ----------
    plane_norm: normal vector of the plane 
    plane_point: a point on the plane 
    ray_direction: direction vector of the ray 
    ray_point: initial position vector of the ray
    
    Output: 
    ---------
    intercept: coordinate of interception between plane and line  
    
    """
    u_dot_n = np.dot(plane_norm, ray_direction)
    
    if abs(u_dot_n) < 1e-6: #can set lower if needed
        return None 

    sep = plane_point - ray_point
    l = np.dot(sep, plane_norm) / u_dot_n
    intercept = ray_point + l * ray_direction 
    return intercept



# Spherical refracting surface

class SphericalRefraction(OpticalElement):
    """
    CLASS USED TO REPRESENT A SPHERICAL SURFACE ON THE Z-AXIS
    
    Attributes: 
    ----------
    z0: the intercept of the surface with the z-axis
    curvature: curvature of the surface
    n1, n2: the refractive indices either side of the surface
    r_aperture: the maximum extent of the surface from the optical axis
    r_aperture should be less than radius 
    
    Methods: 
    ---------
    intercept: calculates intercept of a surface and a ray 
    refract: calculates the direction of the refracted ray from the surface
    propagate_ray: refracts and propagates ray to the intercept point 
    
    """
    def __init__(self, z0, curvature, n1, n2, r_aperture):
        
        self._z0 = z0
        self._curvature = curvature
        self._n1 = n1
        self._n2 = n2
        self._r_aperture = r_aperture 
        
        if self._curvature != 0: #0 curvature would return an infinite radius. 
            self._radius = 1/self._curvature
        else: 
            self._radius = np.inf #Not meant for calculation.
    
        if self._curvature != 0 and self._r_aperture > np.abs(1/self._curvature): 
            raise ValueError('The aperture radius cannot be larger than \
                             the radius of the spherical surface!')
        

       

    def intercept(self, ray): 
        """ 
        FUNCTION THAT CALCULATES INTERCEPT COORDINATE OF A SURFACE AND A RAY 
        VIA VECTOR ANALYSIS.
        
        Parameters: 
        ----------
        
        r: vector from the sphere's centre to the ray's starting point.
        k_hat:  unit vector of the ray.
        p0: initial position of the ray.
        det: determinant of the quadratic eq used to find the intercepts
        centre: vector pointing from the centre of the surface to the 
                ray's intial position
        sol_1: first intercept solution
        sol_2: second intercept solution 
        
        Returns:
        ----------
        
        For a given ray and surface, it returns: 
            The interception coordinate between the ray and the surface
            If the interception point is invalid, the output is 'None'.
        """
        if isinstance(ray, Ray) == False:
             print(isinstance(ray, Ray))
             raise TypeError('the ray input has to be a Ray object.')
        
        centre = np.array([0, 0, 
                           self._z0 + self._radius]) #centre of sphere
        k_hat = ray.k()/np.linalg.norm(ray.k()) #normalise wave vector
        p0 = ray.p() 
        r = p0 - centre

        if self._curvature == 0: #special case: planar surface
            plane_norm = np.array([0, 0, 1])
            plane_point = np.array([0, 0, self._z0])
            intercept = line_plane_intercept(plane_norm, plane_point, 
                                             k_hat, p0)
            return intercept 
        
        
            
        else: #curved surface
        
        #find the possible intersection points
            det = np.dot(r, k_hat)**2 - (np.linalg.norm(r)**2 
                                         - self._radius**2) 
            if det < 0: 
                return None 
            
            else:

                sol1 = -np.dot(r, k_hat) + np.sqrt(det)
                sol2 = -np.dot(r, k_hat) - np.sqrt(det)  
                
                if self._curvature > 0: #convex surface 
                
                    if sol2 > 0:
                        correct_sol = sol2
                        
                    elif sol1 > 0:
                        correct_sol = sol1
                        
                    else:
                        return None
                    
                    I =  correct_sol*k_hat + p0
                    
                    if np.linalg.norm(I[:2]) <= self._r_aperture:
                        return I 
                    
                    else:
                        return None 
                    
                else: #concave surface 
                    
                    if sol1 > 0:
                        correct_sol = sol1
   
                    else:
                        return None
                    
                    I =  correct_sol*k_hat + p0
                    
                    #check intercept is within the aperture 
                    if np.linalg.norm(I[:2]) <= self._r_aperture:
                        return I 
                    
                    else:
                        return None 
                    
              
        
    def propagate_ray(self, ray):
        
        """This method propagates the ray from the starting point to the 
        surface's intercept. 
        
        It then refracts the ray via the refraction method. 

        Finally, it updates the ray's position and direction via 
        the append method. 
        
        If the intercept is invalid, the ray is terminated. 
        
        If total internal reflection occurs inside the surface, the ray 
        is also terminated. 
           """
        intercept = self.intercept(ray)
        if ray._terminate == False: 
            
            if intercept is None:  #check for valid intercept 
                print('No valid intercept. Terminated')
                ray._terminate = True 
            
            
            else:
                
                centre = np.array([0, 0, self._z0 + self._radius])
                p_new = intercept
                
                if self._curvature == 0:
                    n = np.array([0,0 ,-1]) #normal vector of plane in z-axis
                    
                elif self._curvature > 0: #convex surfaces 
                
                    n = (p_new - centre)  #un-normalised 
                    n = n/np.linalg.norm(n) #normalised  
                
                else: #concave surfaces 
                
                    n = -(p_new - centre) 
                    n = n/np.linalg.norm(n)
                    
                k_new = refract(ray.k(), n, self._n1, self._n2)
                
                if k_new is not None: 
                    ray.append(p_new, k_new)

                else: 
                    ray._terminate = True 
            
#Output plane (screen)

class OutputPlane(OpticalElement):
    """
    CLASS USED TO REPRESENT AN OUTPLANE PLANE ON THE Z-AXIS
    
    Attributes: 
    ----------
    z0: position of the output plane on the z-axis
    
    Methods: 
    ---------
    intercept: determines where the ray hits the output plane

    propagate_ray: propagate ray to the output plane without refraction

    """
    def __init__(self, z0):
        self._z0 = z0
        
    def intercept(self, ray):
        
        if ray.k() is None:
            return None
        
        elif ray.k()[2] != 0:
            lamb = (self._z0 - ray.p()[2])/ ray.k()[2] 
            intercept = ray.p() + lamb*ray.k()
            return intercept 
        else: 
            return None
        
    def propagate_ray(self, ray):
        if ray._terminate == False:
            p_new = self.intercept(ray)
            if not p_new is np.ndarray: 
                ray.append(p_new,ray.k()) #ray didn't change direction 

#Beam bundle used for ray tracing investigation 
        
        
        
class Bundle():
        """
        CLASS USED TO REPRESENT A BUNDLE OF RAYS FOR A UNIFORM COLLIMATED BEAM
        ALONG THE Z-AXIS.
        
        Attributes: 
        ----------
        diameter: diameter of the beam
        n: determines # of rays in bundle. n results in n^2 rays 
        (n should be around 10 to avoid long calculations)
        Methods: 
        ---------
        make_bundle: returns n^2 rays uniformly spaced inside the beam.

        """
        def __init__(self, diameter, n):
            self.__diameter = diameter 
            self.__n = n
            if self.__n >= 100: 
                warnings.warn('More than 10^4 rays in beam! \
                             Calculation might take longer!')
            
        def make_bundle(self):
            #r is the radial distance in cylinderical coordinate 
            #phi is the angle phi in cylinderical coordinate
            
            r = np.linspace(-self.__diameter/2, self.__diameter/2, self.__n)
            phi = np.linspace(0,2*np.pi, self.__n)
            rv ,phiv  = np.meshgrid(r,phi) 
            x0 = np.ravel(rv*np.cos(phiv)) #flatten array
            y0 = np.ravel(rv*np.sin(phiv)) #flatten array
            bundle =[]
            direction = [0, 0, 1]
            
            for i in range(len(x0)):
                position = np.array([x0[i], y0[i], 0])
                bundle.append(Ray(position, direction) )
                
            return bundle
        
    

#Aspheric surface

class AsphericRefraction(SphericalRefraction): 
    """
    CLASS USED TO REPRESENT AN ASPHERIC SURFACE ON THE Z-AXIS
    
    ATTRIBUTES 
    ----------
    z0: position of the surface on the z-axis
    
    curvature: curvature of the surface 
    
    A_4, A_6: parameters causing aspheric properties in surface 
    
    n1, n2: refractive index of initial and final medium
    
    r_aperture: aperture radius of surface 
    
    METHODS 
    ----------
    surface_equation: gives equation z(x,y) of aspheric surface
    
    minimise_func: function to be minimised to find intercept between 
                    surface and ray 
                    
    intercept: determines intercept between ray and surface via a minimisation  
               method
               
    normal_vector: calculates the (unit) normal vector of surface pointing
                out of surface. 
    
    propagate_ray: refracts ray at interception point and update ray with   
                    new position and direction vector. 
    
    """
    
    def __init__(self, z0, curvature, A_4, A_6, n1, n2, r_aperture):
        
        self.z0 = z0
        self.curvature = curvature
        self.A_4 = A_4
        self.A_6 = A_6 
        self.n1 = n1
        self.n2 = n2
        self.r_aperture = r_aperture
    

    def surface_equation(self, x, y):

        p = np.sqrt(x**2 + y**2) 
        C = self.curvature 
        A_4 = self. A_4
        A_6 = self.A_6 
        
        term1 = C *p**2/(1+ np.sqrt(1-C*C*p*p))
        term2 = A_4 * p**4
        term3 = A_6 * p**6 
        
        z = self.z0 + term1 + term2 + term3
        
        return z
    
    def minimise_func(self, t,  p0, k):
 
        p_ray = p0 + t*k #point on ray parameterised by t 
        z_value_surface  = self.surface_equation(*p_ray[:2])
        
        return np.sum((z_value_surface - p_ray[2])**2)


    def intercept(self, ray):
        if ray._terminate == False: 
            
            sphere = SphericalRefraction(self.z0, self.curvature, 
                                         self.n1, self.n2, self.r_aperture)
            
            intercept_guess = sphere.intercept(ray)
            
            if type(intercept_guess) is not np.ndarray: 
                ray._terminate = True
            else:
                z_guess = intercept_guess[2]
                
                #intercept occurs when the minimise_func is min. 
                #t_min is when that happens. 
                
                t_min = op.fmin(self.minimise_func, x0=z_guess, 
                                args = ( ray.p(), ray.k() ))
                
                #interection point via line eq
                intersection_point_1 = ray.p() + t_min*ray.k() 
                
                #intersection point via surface equation
                intersection_point_2 = self.surface_equation(
                                            intersection_point_1[0],
                                            intersection_point_1[1])

            #checking that the two points agree + valid
            
                if intersection_point_2 == np.nan: 
                    ray.__terminate = True 
                
                if np.abs(intersection_point_2- intersection_point_1[2])<1e-2:
                    return intersection_point_1
                
                else: 
                    ray._terminate = True 
        
    def normal_vector(self, ray):
        
        if ray._terminate == False: 
            
            C = self.curvature
            A_4 = self.A_4
            A_6 = self.A_6
            
            #(x,y,z) coordinates for the normal to be evaluated at 
            r = self.intercept(ray) 
            
            #Calculate the x y z component of the normal vector
           
            
            p = (r[0]*r[0] + r[1]*r[1])
            term1 = 4* A_4*r[0]*p + 6*A_6*r[0]*p*p 
            
            term2 = C**3 *r[0] *p/( np.sqrt(1-C*C*p)*(1+np.sqrt(1-C*C*p))**2 )
           
            term3 = 2*C*r[0]/(1+np.sqrt(1-C*C*p))
           
         
            x_component = term1 + term2 + term3
            
            y_component = r[1]/r[0] *x_component
            
            z_component = -1 
            
            #putting into array 
            n = np.array( [x_component, y_component, z_component] )
            n = n/np.linalg.norm(n) #unit vector
            
            return n
          
            
        
        
    def propagate_ray(self, ray): 
        
     intercept = self.intercept(ray)
     
     if ray._terminate == False: 
         
         p_new = intercept
         norm = self.normal_vector(ray)
         
         if type(norm) is np.ndarray: 
             
             k_new = refract(ray.k(), self.normal_vector(ray), 
                             self.n1, self.n2)
         
         if type(k_new) is np.ndarray:
             
             ray.append(p_new, k_new)

         else: 
             ray._terminate = True 


#2D lines intersection to determine paraxial focus

def slope(point1, point2):
     return (point2[1]- point1[1])/(point2[0] - point1[0])

def yintercept(point1, point2):
     return point2[1] - slope(point1, point2) * point2[0]

def intercept(ray1, ray2): 
     """
  Function uses slope and yintercept function to determine intercept between
  2 paraxial rays. This allows the paraxial focus to be found.  
    
  Parameters
  ----------
  ray1 : Ray object (praxial ray)
      
  ray2 : Ray object (paraxial ray)

  Returns
  -------
  if rays intercept: z_intercept (float): z coordinate of the ray intercept 
  
  if rays are parallel: returns None
    """
    
    
     #setting up points

     point1_ray1 = np.array( [ray1.p()[2], ray1.p()[0]])
     point2_ray1 = ray1.p() + ray1.k()
     point2_ray1  = np.array( [point2_ray1[2], point2_ray1 [0]])

     point1_ray2 = np.array( [ray2.p()[2], ray2.p()[0]])
     point2_ray2 = ray2.p() + ray2.k()
     point2_ray2  = np.array( [point2_ray2 [2], point2_ray2 [0]])
     
     #calculate slope in z-x plane
     slope1 = slope(point1_ray1, point2_ray1 )
     slope2 = slope(point1_ray2, point2_ray2 )
     
     if slope1 == slope2: 
         print('Parallel rays do not intersect')
         return None
     
     else:   
         xintercept1 = yintercept(point1_ray1, point2_ray1 )
         xintercept2 = yintercept(point1_ray2, point2_ray2 )
     
         z_intercept = (xintercept2 - xintercept1)/(slope1 - slope2)
         # z coordinate of intercept
         
         return z_intercept
 