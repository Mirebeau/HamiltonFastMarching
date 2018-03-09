#!/usr/bin/env python3


"""
    This file demonstrates some EXPERIMENTAL models (not in the base executable):
    - ReedsSheppThreeSpeeds2 and ReedsSheppAdaptive2 model. 
    They both implement the same model: ReedsShepp with (forward, reverse, angular) speed
    specified independently. 
    The implementation of the ReedsSheppAdaptive2 model 
    uses an adaptive angular sampling to guarantee that there will be no 
    sideways motion artifacts.
    - Speed depending only on the angular variable

    Command line (optional) arguments :
    - Path to executable
    - domain size
    - model
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


n = 100 if len(sys.argv)<=2 else int(sys.argv[2])
model = "ReedsSheppThreeSpeeds2" if len(sys.argv)<=3 else sys.argv[3]

# catamaran, wind speed 3
boatSpeed = [0., 0.33, 0.67, 1., 1.33, 1.67, 2., 2., 2., 4.3, 4.7, 5.15, 5.6, \
            5.85, 6.2, 6.35, 6.6, 6.6, 6.7, 6.6, 6.5, 6.6, 6.6, 6.25, 5.9, 5.5, \
            5.1, 4.6, 4.1, 3.6, 3.2, 2.8, 2.4, 2., 2., 2., 2., 2., 2., 2., 2.4, \
            2.8, 3.2, 3.6, 4.1, 4.6, 5.1, 5.5, 5.9, 6.25, 6.6, 6.6, 6.5, 6.6, \
            6.7, 6.6, 6.6, 6.35, 6.2, 5.85, 5.6, 5.15, 4.7, 4.3, 2., 2., 2., \
            1.67, 1.33, 1., 0.67, 0.33]
nTheta = len(boatSpeed)

xi=0.2

if model=="ReedsSheppThreeSpeeds2" or model=="ReedsSheppAdaptive2":
    boatSpeed= [ [x,0,x] for x in boatSpeed]
    """This Reeds-Shepp like model allows to set independently 
    the [forward, reverse, angular] speed."""

input = {
"arrayOrdering": "YXZ_RowMajor", #compatibility with numpy.meshgrid
"model": model,
"dims": np.array([n, n, nTheta]),
"gridScale": 1./n,
"seeds": np.array([[0.5, 0.5, 1.]]),
"exportValues": 1,

"speed":np.roll(boatSpeed,int(round(0.3*nTheta)),axis=0),

"tips":np.array([
                 [0.08, 0.08, 0.0], [0.08, 0.25, 0.0], [0.08, 0.42, 0.0], [0.08, 0.58, 0.0], [0.08, 0.75, 0.0],
                 [0.08, 0.92, 0.0], [0.25, 0.08, 0.0], [0.25, 0.25, 0.0], [0.25, 0.42, 0.0], [0.25, 0.58, 0.0],
                 [0.25, 0.75, 0.0], [0.25, 0.92, 0.0], [0.42, 0.08, 0.0], [0.42, 0.25, 0.0], [0.42, 0.75, 0.0],
                 [0.42, 0.92, 0.0], [0.58, 0.08, 0.0], [0.58, 0.25, 0.0], [0.58, 0.75, 0.0], [0.58, 0.92, 0.0],
                 [0.75, 0.08, 0.0], [0.75, 0.25, 0.0], [0.75, 0.42, 0.0], [0.75, 0.58, 0.0], [0.75, 0.75, 0.0],
                 [0.75, 0.92, 0.0], [0.92, 0.08, 0.0], [0.92, 0.25, 0.0], [0.92, 0.42, 0.0], [0.92, 0.58, 0.0],
                 [0.92, 0.75, 0.0], [0.92, 0.92, 0.0]
                 ]),

"eps":0.1, #does not apply to the ReedsSheppForwardAdaptive2 model
"xi":xi,
"uniformlySampledSpeed":nTheta,
}

# Get the executable name and path
FileHFM_executable = "FileHFM_"+model
if len(sys.argv) >=2:   FileHFM_binary_dir = sys.argv[1];

#Write parameter file to disk, execute program, import results
result = FileIO.WriteCallRead(input,
                              FileHFM_executable, binary_dir=FileHFM_binary_dir,
                              logFile=None) # Display output directly

geodesics = np.vsplit(result['geodesicPoints'],result['geodesicLengths'].astype(int).cumsum()[:-1])

with PdfPages("Output/Boat_"+model+"_results.pdf") as pdf:
    plt.figure()
    for geo in geodesics:
        plt.plot(geo[:,0],geo[:,1],'r')
    plt.title('Some geodesics')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
 
