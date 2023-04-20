#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Required input modules, global variable declarations and
helper functions for SIF_final.py
"""
# ------------ required modeules --------------------------------- #
import math
import numpy as np
import matplotlib.pyplot as plt
from numpy import random
#from numpy.core.umath_tests import inner1d
# from matplotlib.tri import Triangulation
import shapely.geometry as geometry
from shapely.geometry import LineString, MultiLineString
from scipy.spatial import Delaunay
# from matplotlib.collections import LineCollection
from matplotlib.path import Path
import pdb ### FOR DEBUGGING ###
# import lmp_class
from psifa.api import psifa #developed by Ferreira, Ehrhardt and Santos (2017)



# ------------ global variables --------------------------------- #

pi = math.pi
R_0 = 7                  # max radius for triangulization routine 
poisson =  0.25          # material specific
E = 70                   # material specific Youngs modulus
kolosov = 3-4*poisson    # material specific      
shear_modulus = 0.25     # material specific  
random_or_all = 'a'      # chose points randomly from the annulus or us all?
nmax = 2                  # number of terms used in the series expansion
nrpoints = 500           # if random (above), how many points do you choose? 
R1 = 50  #R1 = 20;       # min size for annulus:  analyse the displacment field in an annulus defined by R1 and R2
R2 = 112#R2 = 120;       # max size for annulus:  analyse the displacment field in an annulus defined by R1 and R2... must be <= 150
boxsize = 4*R2             # defines a region of interest.  updated cracktip is assumed to fit within a square of this size
tmin = -pi                # minimum theta angles:
tmax = pi                 # maximum theta angles:
N = 10000                 # number of points to use in the initial toy model configuration  
# N = np.linspace(5000,505000,51) #to compare runtimes for different sizes if the system/geometry 
# tip_guess = [10,10];     # Patternsearch initial tip guess
# nnodes = 9801;           #Nodes imported from FEM_result Text File

# ------------ helper functions --------------------------------- #

def cart2pol(x, y):
    ###################################################################
    # # # Helper function # # #
    # Convert cartesian coordinates to polar coordinates
    ###################################################################
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return [theta, r]

def pol2cart(theta, r):
    ###################################################################
    # # # Helper function # # #
    # Convert polar coordinates to cartesian coordinates
    ###################################################################
    
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return [x, y]

def iscolumn(x):
    #################################################################
    # HELPER FUNCTION
    #
    # This function replicates functionality of 
    # Matlab's "iscolumn" that "returns logical 1 (true) if size(V) 
    # returns [n 1] with a nonnegative integer value n, and 
    # logical 0 (false) otherwise."
    ##################################################################
    sizeX = np.shape(x);
    try:
        check = sizeX[0] >= 1 and sizeX[1] == 1;
        # print (check);
        # print (sizeX);
    except IndexError:
        return False;
    if (check):
        return True;
    else:
        return False;
        
def SVD_solution(matr, vect):
    
    
    U,S,V = np.linalg.svd(matr,0)
    
    x = np.linalg.inv(U @ S @ V.T) @ vect
    
    return x     
   
def randperm(n,k):
############################################################################
# # # Helper function # # #
#
# randperm(n,k) returns a row vector containing K unique integers
# selected randomly from 1:N.
############################################################################
    r = range(0, n);

    lst = [i for i in r];
   
    random.shuffle(lst)
    ids = range(0,k)
    rp = [lst[i] for i in ids]

    return rp;

def inpolygon(x,y,xv,yv):
    ############################################################################
    # HELPER Function
    # Implementtaion of Matlab's inpolygon(X,Y,XV,YV)
    ############################################################################

    l = len(xv)
    vertices = [(xv[i],yv[i]) for i in range(0,l)]
    path = Path(vertices)
    testPoints = np.hstack([x.reshape(x.size, -1,order = 'F'), y.reshape(y.size,-1,order = 'F')])
    _in = path.contains_points(testPoints)
    _in_on = path.contains_points(testPoints, radius = 1e-01)
    
    _on = _in ^ _in_on
    
    return _in_on, _on


def determine_points_of_interest_and_displacement(lmp,lmp_reference, cracktip2,config):
    ############################################################################
    #
    # Function will calculate a displacement field based on two different MD snapshots. 
    # lmp is a current snapshot and lmp_reference is a snapshot earlier
    # in time. The function identifies atom ids, dx, dy, r and theta values for 
    # a desired number of points (atoms) within a defined annulus size. 
    #   
    # function used in MD analysis only
    #    
    ############################################################################
    
    pi = math.pi;
    ## no point of calculating r here as its recalculated in cart2pol
    #r = (cracktip2[0]-lmp.coordinate[:,0])**2 + (cracktip2[1]-lmp.coordinate[:,1])**2;
    
    ### TODO: is calculating theta here redundant?
    
    if config:
        [theta,r] = cart2pol(lmp.coordinate[:,0]-cracktip2[0], lmp.coordinate[:,1]-cracktip2[1]);
        theta = theta + pi/2;
        theta[theta > pi] =  theta[theta > pi] - 2*pi;                    
        Ax = [lmp.coordinate[:,0], lmp.coordinate[:,1]];
        Bx = [lmp_reference.coordinate[:,0], lmp_reference.coordinate[:,1]];
    else:
        ids_SiO = np.where((lmp.type == 1) | (lmp.type == 2))[0]
        ids_SiO_reference = np.where((lmp_reference.type == 1) | (lmp_reference.type == 2))[0]
        
        [theta,r] = cart2pol(lmp.coordinate[ids_SiO,0]-cracktip2[0], lmp.coordinate[ids_SiO,1]-cracktip2[1])
        theta = theta + pi/2;
        theta[theta > pi] =  theta[theta > pi] - 2*pi;                    
        Ax = [lmp.coordinate[ids_SiO,0], lmp.coordinate[ids_SiO,1]];
        Bx = [lmp_reference.coordinate[ids_SiO_reference,0], lmp_reference.coordinate[ids_SiO_reference,1]]
    
    DX = np.array([Ax])- np.array([Bx]); 
    
    return DX, r, theta;

def gen_box_poly(boxsize, x0, y0):
    #####################################################################
    # This function builds a polygon around the general 
    # know region of interest 
    #####################################################################
    xmin = x0-boxsize/2;
    xmax = x0+boxsize/2;
    ymin = y0-boxsize/2;
    ymax = y0+boxsize/2;
    boxx = [xmin, xmin, xmax, xmax, xmin, np.nan];
    boxy = [ymin, ymax, ymax, ymin, ymin, np.nan];

    return boxx,boxy;


def find_boundary(X):
    #####################################################################
    # This function returns the boundary B indicies of points that mark the 
    # boundary of the (x,y) locations of a 2D set of points 
    #####################################################################
    plotBoundary = 0;
    lengthX = max(np.shape(X))
    points = np.array( [ [ X[0,i] , X[1,i]] for i in range(0,lengthX) ] )

    edges = alpha_shape(points, alpha=R_0, only_outer=True)
    
    # fig, ax = plot_polygon(concaveHull)
    if(plotBoundary):
        fig, ax = plt.subplots()
        for i, j in edges:
            ax.scatter(points[[i, j], 0], points[[i, j], 1], marker = '.', color = "blue")
            
        ax.set_xlabel("X Direction")
        ax.set_ylabel("Y Direction")
        ax.axis('square')
        ax.set_title('Boundary Points using Delaunay Triangulation ')
        # plt.plot(points,'o')
        # ax.scatter(edges)
        
        plt.show()
    
    # edges is a set (i.e unique values only)
    return edges;

# modified by SG from :
# https://stackoverflow.com/questions/23073170/calculate-bounding-polygon-of-alpha-shape-from-the-delaunay-triangulation
def alpha_shape(points, alpha, only_outer=True):
    """
    Compute the alpha shape (concave hull) of a set of points.
    :param points: np.array of shape (n,2) points.
    :param alpha: alpha value.
    :param only_outer: boolean value to specify if we keep only the outer border
    or also inner edges.
    :return: set of (i,j) pairs representing edges of the alpha-shape. (i,j) are
    the indices in the points array.
    """
    assert points.shape[0] > 3, "Need at least four points"

    def add_edge(edges, edge_points, coords,  i, j):
        """
        Add a line between the i-th and j-th points,
        if not in the list already
        """
        if (i, j) in edges or (j, i) in edges:
            # already added
            assert (j, i) in edges, "Cannot go twice over same directed edge"
            if only_outer:
                # if both neighboring triangles are in shape, it is not a boundary edge
                edges.remove((j, i))
            return
        edges.add((i, j))
        edge_points.append(coords[ [i, j] ])
    
    
    #coords = np.array([point[0]for point in points])
    tri = Delaunay(points)
    edges = set()
    edge_points = []
    # Loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in tri.simplices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        # Computing radius of triangle circumcircle
        # www.mathalino.com/reviewer/derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle
        # Lengths of sides of triangle
        a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        # Semiperimeter of triangle
        s = (a + b + c) / 2.0
        # Area of triangle by Heron's formula
        area = np.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)
        # radius filter
        if circum_r < alpha:
            add_edge(edges,edge_points, points,  ia, ib)
            add_edge(edges,edge_points, points,  ib, ic)
            add_edge(edges,edge_points, points,  ic, ia)
    
    m = geometry.MultiLineString(edge_points)
    return edges 
        
 