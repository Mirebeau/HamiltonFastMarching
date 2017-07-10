% Copyright Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay

% This file demonstrate anisotropic fast marching with the "Half-Disk" model.
% Since it is not among the "standard" models, you need to build the executable MatlabHFM_RiemannExtra, using the script CompileMexHFM.m

clear input;
nx=81;ny=83; %Taking different dimensions to make sure axes are not inverted...
nGeo=8;

% 'Trivial' means that we work on plain R^d, without a bundle-like component.
input.model = 'HalfDisk2';

% We do use Matlab's convention of ordering the axes as y, x, z.
% Internally, the first two array coordinates are transposed.
% This feature can be turned off by setting input.transposeFirstTwoCoordinates=False
[x,y]=meshgrid(1:nx,1:ny);

%   dualMetric [Vx,Vy,r] for the hamiltonian V o V + r^2 Vp o Vp, where V=(Vx,Vy), and Vp is the orthogonal vector.
% In other words, the front goes in the direction of V at speed |V| (the norm of V), and in the sideways direction of +Vp and -Vp at speed r |V|.
% Motion in the direction of -V is in principle forbidden (in practice, it still happens at speed eps r |V|).
input.dualMetric = zeros([3,ny,nx]);
input.dualMetric(1,:,:)=1;      % Vx
input.dualMetric(2,:,:)=0.5;    % Vy
input.dualMetric(3,:,:)=0.5;    % r

% Relaxation parameters are:
% - a scalar eps, for how badly is backwards motion penalized.
% - a scalar epsForward, for how accurately forward speed is implemented.
% As usual, small parameters yield large stencils. We leave default values here.

% The three dimensional model HalfDisk3 is completely similar, except that the dualMetric format reads [Vx,Vy,Vz,r].


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

output=MatlabHFM_RiemannExtra(input);

clf;
imagesc(output.values)
geodesics = mat2cell(output.geodesicPoints,2,output.geodesicLengths);
for i=1:size(geodesics,2) %geodesics joining the "tips"
    rescaledGeodesic=RescaledCoords(geodesics{i},input.origin,[input.gridScale;input.gridScale]);
    line(rescaledGeodesic(1,:),rescaledGeodesic(2,:));
end;
pause;