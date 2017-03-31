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


Program description:
This code implements a Fast-Marching solver with adaptive stencils. It computes distance maps and shortest paths with respect to a variety of metrics, on a domain discretized on a cartesian grid. Typical uses include motion planning and image segmentation. The code features:
- riemannian metrics, in dimension 2 and 3, and curvature penalizing metrics such as the Reeds-Shepp, Euler Elastica or Dubins models.
- second order accuracy (optional), various stopping criteria, propagation of states, forward and backward differentiation.
- interfaces to  Mathematica(R) and Python(R) using files, Matlab(R) using mex, see ExampleFiles directory.

Installation instructions:
-You also need to download the github repository github.com/Mirebeau/JMM_CPPLibs (header only).
-The code is written in C++11. Use CMAKE to build, or Matlab(R) mex. You may need need to enter the path to to JMM_CPPLibs.

Fair use:
If you use this program for an academic or commercial project, then please cite at least one of the following papers.

- (On the Reed-Shepp models) : 
Remco Duits, Stephan P.L. Meesters, Jean-Marie Mirebeau, Jorg M. Portegies, Optimal Paths for Variants of the 2D and 3D Reeds-Shepp Car with Applications in Image Analysis, Preprint available on arXiv, 2017

- (On Riemannian metrics) : 
Jean-Marie Mirebeau, Anisotropic fast-marching on cartesian grids using Voronoi's first reduction of quadratic forms, in preparation, 2017

- (On Euler Elastica and Dubins paths) : 
Jean-Marie Mirebeau, Fast Marching methods for Curvature Penalized Shortest Paths, in preparation, 2017