#!/usr/bin/env python3


"""
    This file lets you compute geodesics in a square domain with a wall.
    Command line (optional) arguments :
    - Path to executable
    - domain size
    - model (ReedsShepp2, ReedsSheppForward2, Elastica2<5>, Dubins2)
"""

import math
import numbers;
import numpy as np;
import os;
import sys
from operator import mul
from functools import reduce
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import FileIO


size = 61 if len(sys.argv)<=2 else int(sys.argv[2])
model = "Elastica2<5>" if len(sys.argv)<=3 else sys.argv[3]
#model = 'Elastica2<5>' #'Dubins2' #'ReedsShepp2'


def SquareTest(n,model,withWall=False):
    input={
        'eps':0.1,
        # space0,space1,angle
        'dims':np.array([n+1,n,60]),
        # size of a pixel (only for physical dimensions)
        'gridScale':1./n,
        'model':model,
        # Speed(s) multiplier
        # constant
        # 'multiplier':1,
        # array (size of the physical vector)
        # array (size of angular space)
        # array of size the space
        #'multiplier':np.array([n,n,60],1), 
        'speed':1,

        # typical radius of curvature
        'xi': 0.1,
        
        # Starting point(s) for the front propagation.
        # pos_x, pos_y, angle_theta
        'seeds':np.array([[0.5,0.5,1.]]),
        # Alternatively (or in addition), omni-directional starting points can be specified. # pos_x, pos_y
        # 'seeds_Unoriented':np.array([[0.5,0.5]]),
        
        # Starting points for geodesic backtracking
        # pos_x, pos_y, angle_theta
        'tips':np.array([
            [0.08, 0.08, 0.0], [0.08, 0.25, 0.0], [0.08, 0.42, 0.0], [0.08, 0.58, 0.0], [0.08, 0.75, 0.0],
            [0.08, 0.92, 0.0], [0.25, 0.08, 0.0], [0.25, 0.25, 0.0], [0.25, 0.42, 0.0], [0.25, 0.58, 0.0],
            [0.25, 0.75, 0.0], [0.25, 0.92, 0.0], [0.42, 0.08, 0.0], [0.42, 0.25, 0.0], [0.42, 0.75, 0.0],
            [0.42, 0.92, 0.0], [0.58, 0.08, 0.0], [0.58, 0.25, 0.0], [0.58, 0.75, 0.0], [0.58, 0.92, 0.0],
            [0.75, 0.08, 0.0], [0.75, 0.25, 0.0], [0.75, 0.42, 0.0], [0.75, 0.58, 0.0], [0.75, 0.75, 0.0],
            [0.75, 0.92, 0.0], [0.92, 0.08, 0.0], [0.92, 0.25, 0.0], [0.92, 0.42, 0.0], [0.92, 0.58, 0.0],
            [0.92, 0.75, 0.0], [0.92, 0.92, 0.0]
        ]),
        # Alternatively (or in addition), omni-directional tips can be specified. The minimal orientation is automatically selected. # pos_x, pos_y
        'tips_Unoriented':
        np.array([
            [0.08, 0.08], [0.08, 0.25], [0.08, 0.42], [0.08, 0.58], [0.08, 0.75],
            [0.08, 0.92], [0.25, 0.08], [0.25, 0.25], [0.25, 0.42], [0.25, 0.58],
            [0.25, 0.75], [0.25, 0.92], [0.42, 0.08], [0.42, 0.25], [0.42, 0.75],
            [0.42, 0.92], [0.58, 0.08], [0.58, 0.25], [0.58, 0.75], [0.58, 0.92],
            [0.75, 0.08], [0.75, 0.25], [0.75, 0.42], [0.75, 0.58], [0.75, 0.75],
            [0.75, 0.92], [0.92, 0.08], [0.92, 0.25], [0.92, 0.42], [0.92, 0.58],
            [0.92, 0.75], [0.92, 0.92]
        ]),
        'verbosity':2,
    }
    if withWall:
        input['walls']=np.array([[i==math.floor(n/3) and j<2*n/3 for j in range(n)] for i in range(n+1)])
    # This defines a wall depending on the physical coordinates. Alternatively one may define a wall depending on both physical and bundle coordinates.
    return input

withWall=True;

# Set to the directory containing ExternLFM executables FileLFM2D and FileLFM3D
path=os.getcwd()
if len(sys.argv) >=2:
    os.chdir(sys.argv[1])

#Write parameter file to disk, execute program, import results
FileIO.RulesToRaw(SquareTest(size,model,withWall),"input")
# call(path+'/FileHFM_Curvature2 | tee log.txt')
call('./FileHFM_AllBase') # Also in FileHFM_Curvature2
result = FileIO.RawToRules()

# print the result to see additional output fields
#print(result)

# Restore working directory
os.chdir(path)
geodesics = np.vsplit(result['geodesicPoints'],result['geodesicLengths'].astype(int).cumsum()[:-1])
geodesicsU = np.vsplit(result['geodesicPoints_Unoriented'],result['geodesicLengths_Unoriented'].astype(int).cumsum()[:-1])

with PdfPages("Output/Curvature2_"+model+'_results.pdf') as pdf:
    plt.figure()
    #plt.plot(range(7), [3, 1, 4, 1, 5, 9, 2], 'r-o')
    #plt.plot([1,2,4,3],[2,5,4,2],'r')
    for geo in geodesics:
        plt.plot(geo[:,0],geo[:,1],'r')
    if withWall:
        plt.plot([1/3,1/3],[0,2/3],'b')
    plt.title('Some geodesics')
    pdf.savefig()  # saves the current figure into a pdf page

    plt.figure()
    for geo in geodesicsU:
        plt.plot(geo[:,0],geo[:,1],'r')
    if withWall:
        plt.plot([1/3,1/3],[0,2/3],'b')
    plt.title('Some geodesics, with unoriented tips')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()


 
