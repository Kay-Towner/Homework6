#Finished by: Kay Towner

import numpy as np
import math
import scipy.integrate as integrate
import numpy.ma as ma

#Problem 2.) Flat Earth

def grav_accel(p1, p2, m):
    """ p1 = point where the mass element is
        p2 = point you are interested in
        m  = mass
        returns a vector of the gravitational accleration"""
    G = 6.6743e-11
    diff = p1 - p2
    r = np.sqrt(np.sum(diff))
    rhat = diff / r  
    return -1*G*m/r**2*rhat


def point_in_flatearth(x, y, z, radius=None):
    """Function for the Mask Condition."""
    #FILL THIS IN WITH A VALID MASK CONDITION
    #THAT RETURNS TRUE WHEN X,Y,Z IN THE SPHERE
    if (x - center_x)**2 <= (radius**2):
        return True
    if (y - center_y)**2 <= (radius**2):
        return True
    if (z - center_z)**2 <= (radius**2):
        return True
    else:
        return False


if __name__ == "__main__":
    km = 1000 #1 km = 1000 meters
    rho = 5514 #kg/m^3, density of flat Earth
    r_flat_earth = 20037*km #radius of flat Earth
    h = 200.0*km #relatively coarse step size
    dV = (1/3)*np.pi*(r_earth)**3#volume of sphere
    #set grid size same in x,y,z
    dx, dy, dz = h, h, h
    #x, y, z define boundaries of grid, here 7000 km
    x = np.arange(-7000*km, 7000*km, dx)
    len_x = x.shape[0]
    y = x.copy()
    z = y.copy()
    
    #define points on the north pole, south pole, and equator
    point_northpole = np.array([0, 0, 20037*km])
    point_equator   = np.array([20037*km, 0, 0])
    point_southpole = np.array([0, 0, -20037*km])

    ##Sample North Pole calc
    ##(I just copy and pasted the code for each of the points:)
    grav_vec_northpole = np.array([0, 0, 20037*km])
    for idx, xx in enumerate(x):
        #this is a trick to tell how long it will take
        print(idx, " of ", len_x, "x steps.")
        for yy in y:
            for zz in z:
                if point_in_sphere(xx,yy,zz):#FIX THIS, MAKE THIS A VALID FUNCTION
                    m = rho * dV
                    point = np.array([xx, yy, zz])
                    grav_vec_northpole += grav_accel(point, point_northpole, m)                 
    print("The gravity vector at the north pole is...", grav_vec_northpole)
    print("Should be something like [0,0,-9.8] m/s^2")
   
        ##Sample Equator calc
    grav_vec_equatior = np.array([20037*km, 0, 0])
    for idx, xx in enumerate(x):
        #this is a trick to tell how long it will take
        print(idx, " of ", len_x, "x steps.")
        for yy in y:
            for zz in z:
                if point_in_sphere(xx,yy,zz):#FIX THIS, MAKE THIS A VALID FUNCTION
                    m = rho * dV
                    point = np.array([xx, yy, zz])
                    grav_vec_equator += grav_accel(point, point_equator, m) 
         print("The gravity vector at the equator is...", grav_vec_equator))

        ##Sample South Pole calc
    grav_vec_southpole = np.array([0, 0, -20037*km])
    for idx, xx in enumerate(x):
        #this is a trick to tell how long it will take
        print(idx, " of ", len_x, "x steps.")
        for yy in y:
            for zz in z:
                if point_in_sphere(xx,yy,zz):#FIX THIS, MAKE THIS A VALID FUNCTION
                    m = rho * dV
                    point = np.array([xx, yy, zz])
                    grav_vec_southpole += grav_accel(point, point_southpole, m)
        print("The gravity vector at the south pole is...", grav_vec_southpole)




                    
