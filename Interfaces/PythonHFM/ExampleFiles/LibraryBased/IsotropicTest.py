#!/usr/bin/env python

import numpy as np

from HFMpy import HFM_Isotropic2

hfm = HFM_Isotropic2.HFMIO()

from HFMpy import HFM_Riemann2
hfm2 = HFM_Riemann2.HFMIO()

# Demonstrating basic input output

hfm.set_string("Hi","There")
print("Hi:",hfm.get_string("Hi"))

hfm.set_scalar("booh",2)
print("booh:",hfm.get_scalar("booh"))

#arr=np.array([[1.,2.,3.],[4.,5.,6.]]) #also works
arr=np.array([[1,2,3],[4,5,6]]).astype(float)
print("arr:",arr)
hfm.set_array("hum",arr)
arr2=hfm.get_array("hum")
print("arr2:",arr2)

# Setting up a shortest path problem

n=5
hfm.set_array("dims",np.array([2*n,n]).astype(float))

x,y=np.meshgrid(np.linspace(0,2,2*n),np.linspace(0,1,n))
hfm.set_array("speed",np.exp(-(x**2+y**2)))

hfm.set_scalar("gridScale",1./n)
hfm.set_array("seeds",np.array([[0.5,0.5]]))
hfm.set_scalar("exportValues",1)
hfm.set_scalar("verbosity",2)
hfm.run()
values=hfm.get_array("values")

print("values.shape",values.shape)
print("values:",values)
print(hfm.computed_keys())
