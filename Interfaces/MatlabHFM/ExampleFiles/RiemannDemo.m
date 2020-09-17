% Copyright Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay

% This file demonstrate riemannian anisotropic fast marching.

% If you have not compiled this model yet, please run the CompileMexHFM script and then
% compileModelsHFM(binary_Dir,{'Riemann2'})

clear input;
nx=81;ny=83; %Taking different dimensions to make sure axes are not inverted...
nGeo=8;

% We do use Matlab's convention of ordering the axes as y, x, z.
% Use the input.arrayOrdering field to change this behavior
[x,y]=meshgrid(1:nx,1:ny);


%    Two dimensional :  xx, xy, yy
%    Three dimensional  xx, xy, yy, xz, yz, zz
input.dualMetric = zeros([3,ny,nx]);
input.dualMetric(1,:,:)=1;
input.dualMetric(2,:,:)=0.5*((x>nx/2)-(x<nx/2));
input.dualMetric(3,:,:)=1+(y>0.75*ny);
% (In contrast with this example, convergence analysis was only done for continuous metrics.)

input.origin=[0;0]; % grid position parameters
input.gridScale=1/nx;
input.dims = [nx;ny];
%        input.stopWhenTipsReached=1;

input.seeds=[0.5;0.5]; % Where the geodesics start [x;y]. You can set multiple seeds, as for the tips below.

epsGeo = 1/nGeo;
arrGeo = (epsGeo/2):epsGeo:(1-epsGeo/2);
[x,y]=meshgrid(arrGeo,arrGeo);
x=reshape(x,[numel(x),1]); y=reshape(y,[numel(y),1]);
input.tips = [x,y]'; % where the geodesics end. This will select the angle t for which the geodesic is the shortest.
input.exportValues=1; % distance table, of size [n,n] (minimum of the previous one over directions)

output=MatlabHFM_Riemann2(input);

clf;
imagesc(output.values)
geodesics = mat2cell(output.geodesicPoints,2,output.geodesicLengths);
for i=1:size(geodesics,2) %geodesics joining the "tips"
    rescaledGeodesic=RescaledCoords(geodesics{i},input.origin,[input.gridScale;input.gridScale]);
    line(rescaledGeodesic(1,:),rescaledGeodesic(2,:));
end