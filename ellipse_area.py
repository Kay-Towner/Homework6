#by Kay Towner
#Problem1: Area of an Ellipse:

import numpy as np
import math
import matplotlib.pyplot as py
from scipy import integrate
import numpy.ma as ma


def ellipse_equation(x, a=0, b=0):
    """Function of the ellipse equation:
    (x**2/a^2) + (y^2/b^2) = 1
    a = 2
    b = 4 """ 
    return b**2*np.sqrt(1-x**2/a**2)

def area(a, b):
    """Function for the ellipse area formula"""
    return np.pi * a * b
    
#1.)1D

if __name__ == "__main__":
    a = 2
    b = 4
    n = 100
    print("This is the area of the ellipse:",area(a,b))

    xrange = np.linspace(0,a,n) 
    #1D Array:
    simpson = integrate.simps(ellipse_equation(xrange,a,b), xrange)
    print("The area of the ellipse by the Simpsons rule:", simpson)

#2.)2D
    #Mask section:
    yrange = np.linspace(0,b,n)
    mask = np.meshgrid(xrange, yrange)
    mask = ellipse_equation(xrange, a, b) <= 1
    print("Plot Evaluated:", mask)
    









        
