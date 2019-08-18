    HamiltonFastMarching - A Fast-Marching solver with adaptive stencils.

    Copyright (C) 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay.  
    jm 'dot' mirebeau 'at' gmail 'dot' com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Hamiltonian Fast Marching
## A Fast-Marching solver with adaptive stencils.


## Latest version

The latest version of this program is available on the Github repository [HamiltonFastMarching](https://github.com/Mirebeau/HamiltonFastMarching). See also the illustrating [notebooks](http://nbviewer.jupyter.org/urls/rawgithub.com/Mirebeau/AdaptiveGridDiscretizations/master/Summary.ipynb).

The present version is part of a submission to the [Image Processing On Line](http://www.ipol.im) journal. It was validated as reproducible research, and supplemented with an online [demo](http://ipol-geometry.loria.fr/~kerautre/ipol_demo/DemoIPOL_LORIA-Geom_HFM/).

## Program description
This code implements a Fast-Marching solver with adaptive stencils. It computes distance maps and shortest paths with respect to a variety of metrics, on a domain discretized on a cartesian grid. Typical uses include motion planning and image segmentation. The code features:
- Riemannian metrics, in dimension 2 and 3, and curvature penalizing metrics such as the Reeds-Shepp, Euler Elastica or Dubins models.
- second order accuracy (optional), various stopping criteria, propagation of states, forward and backward differentiation.
- interfaces to  Mathematica(R) and Python(R) using files, Matlab(R) using mex, see ExampleFiles directory.

## Installation instructions
- You also need to download the Github repository github.com/Mirebeau/JMM_CPPLibs. When compiling, you may need need to enter the path to this (header only) library.
- The code is written in C++17.
- It is meant to be called from one of the following scripting languages.
  * Python(R). Please compile the code located in Interfaces/FileHFM using CMake. Then consult the examples located in Interfaces/PythonHFM/ExampleFiles/FileBased as well as the illustrative notebooks.
  * Mathematica(R). Please compile the code located in Interfaces/FileHFM using CMAKE. Then consult the examples located in Interfaces/MathematicaHFM/ExampleFiles/FileBased as well as the illustrative notebooks.
  * Matlab(R). Please go to directory Interfaces/MatlabHFM/ExampleFiles in Matlab(R). Then use the script CompileMexHFM to build the executables, and consult the provided examples.

### TL-DR. Simple and direct compilation instructions for Python usage (tested on macOs)
* Extract the archives *HamiltonFastMarching* and *JMM_CPPLibs* in a common directory.
* Open a terminal, cd to the directory containing the extracted archives, and type the following commands.
 ```console
 mkdir bin
 cd bin
 cmake ../HamiltonFastMarching-master/Interfaces/FileHFM
 make
 ```
* Have a look at the [notebooks](http://nbviewer.jupyter.org/urls/rawgithub.com/Mirebeau/AdaptiveGridDiscretizations/master/Summary.ipynb).

## Fair use
If you use this program for an academic or commercial project, then please cite at least one of the following papers.

- (On the Reed-Shepp models) :
Remco Duits, Stephan P.L. Meesters, Jean-Marie Mirebeau, Jorg M. Portegies, Optimal Paths for Variants of the 2D and 3D Reeds-Shepp Car with Applications in Image Analysis, Preprint available on arXiv, 2017

- (On Euler Elastica and Dubins paths) :
Jean-Marie Mirebeau, Fast Marching methods for Curvature Penalized Shortest Paths, Journal of Mathematical Imaging and Vision, 2017

- (On Riemannian metrics) :
Jean-Marie Mirebeau, Anisotropic fast-marching on cartesian grids using Voronoi's first reduction of quadratic forms, Preprint available on HAL, 2018


## Reproducible research
An earlier version of this code was part of a reproducible research publication:
Jean-Marie Mirebeau, Jorg M. Portegies, Hamiltonian Fast Marching: A numerical solver for anisotropic and non-holonomic eikonal PDEs, 2019, accepted for publication.
