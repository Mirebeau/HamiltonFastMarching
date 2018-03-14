% Copyright Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay

% This script demonstrates the propagation of a state in the FM algorithm.
%- In line comutation of geodesic euclidean length.
%- stopping criterion based on geodesic euclidean length
%- voronoi diagram computation

% If you have not compiled this model yet, please run the CompileMexHFM script and then
% compileModelsHFM(binary_Dir,JMM_CPPLibs_Dir,{'Isotropic2'})

n=20;h=1/n;
[x,y]=meshgrid(0:1/n:2,0:1/n:1); 

params.dims=size(x')';
params.speed=1+(x>1);
params.seeds = [0.2,0.8;0.5,0.7];
params.exportValues=1;
params.gridScale=h;

params.euclideanScale = h*[1;1];
params.stopAtEuclideanLength=1;

params.seedFlags=[1;2];

out=MatlabHFM_Isotropic2(params);

imagesc(out.values); title('Seed distance value');
pause;
imagesc(out.voronoiFlags); title('Voronoi diagram segmentation');
pause;
imagesc(out.euclideanLengths); title('Geodesic euclidean length');