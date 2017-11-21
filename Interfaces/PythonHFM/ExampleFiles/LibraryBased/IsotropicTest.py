#!/usr/bin/env python

import numpy as np
import sys

"""
    This file demonstrates basic input output with the PythonHFM_AllBase library
    The library must be visible, or the path to its directory provided as first argument to this script.
"""


if len(sys.argv) >=2:
    sys.path.insert(0,sys.argv[1])

import PythonHFM_AllBase
hfm = PythonHFM_AllBase.HFMIO()

# Demonstrating basic input output

hfm.SetString("Hi","There")
print("Hi:",hfm.GetString("Hi"))

hfm.SetScalar("booh",2)
print("booh:",hfm.GetScalar("booh"))

#arr=np.array([[1.,2.,3.],[4.,5.,6.]]) #also works
arr=np.array([[1,2,3],[4,5,6]]).astype(float)
print("arr:",arr)
hfm.SetArray("hum",arr)
arr2=hfm.GetArray("hum")
print("arr2:",arr2)

# Setting up a shortest path problem

n=5
hfm.SetArray("dims",np.array([2*n,n]).astype(float))

x,y=np.meshgrid(np.linspace(0,2,2*n),np.linspace(0,1,n))
hfm.SetArray("speed",np.exp(-(x**2+y**2)))

hfm.SetString("model","IsotropicBox2<Boundary::Closed>")
hfm.SetScalar("gridScale",1./n)
hfm.SetArray("seeds",np.array([[0.5,0.5]]))
hfm.SetScalar("exportValues",1)
hfm.SetScalar("verbosity",2)
hfm.Run()
values=hfm.GetArray("values")

print("values.shape",values.shape)
print("values:",values)
print(hfm.GetComputedKeys())
