// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

# Tensor decomposition by Voronoi's first reduction

This c++ program is used to decompose positive definite tensors using Voronoi's first reduction method.
It relies on the header only library [JMM_CPPLibs](https://github.com/Mirebeau/JMM_CPPLibs)
Input and output are through a file. (Direct interfaces may be added in the future.)

These techniques are at the foundation of adaptive schemes for anisotropic PDEs on cartesian grids, presented in
[AGD_Python_Notebooks](https://github.com/Mirebeau/AGD_Python_Notebooks)
