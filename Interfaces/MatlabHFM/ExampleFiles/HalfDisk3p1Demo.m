% Copyright Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay

% This file demonstrates anisotropic fast marching with the "Half-Disk" model, in dimension 3+1.
% Three physical dimensions, and one extra isotropic dimension.
% Since it is not among the "standard" models, you need to build the executable MatlabHFM_RiemannExtra, using the script CompileMexHFM.m

clear input;

nx=20;
input.model = 'HalfDisk3p1';
input.origin=[0;0;0;0]; % grid position parameters
input.gridScale=1/nx;
input.dims = [nx;nx+1;nx+2;nx+3];
input.sndOrder=1;
% TODO : distinct grid scale for the physical and bundle coordinates.

nGeo=8;

% Coordinates can be ordered in : Transposed manner y,x,z,t, or Default manner x,y,z,t.
% When calling from the Matlab interface, the default setting is ... Transposed.
input.arrayOrdering='Transposed'; input.dualMetric = zeros([3+1+1,nx+1,nx,nx+2,nx+3]); % Actual default from Matlab interface.
%input.arrayOrdering='Default'; input.dualMetric = zeros([3+1+1,input.dims']);

%   dualMetric [Vx,Vy,r] for the hamiltonian V o V + r^2 Vp o Vp, where V=(Vx,Vy), and Vp is the orthogonal vector.
% In other words, the front goes in the direction of V at speed |V| (the norm of V), and in the sideways direction of +Vp and -Vp at speed r |V|.
% Motion in the direction of -V is in principle forbidden (in practice, it still happens at speed eps r |V|).
% last component is the bundle speed.
input.dualMetric(1,:,:,:,:)=1;      % Vx
input.dualMetric(2,:,:,:,:)=0;      % Vy
input.dualMetric(3,:,:,:,:)=0;      % Vz
input.dualMetric(4,:,:,:,:)=1;      % r
input.dualMetric(5,:,:,:,:)=1.;     % bundle speed

input.seeds=[0.5;0.5;0.5;0.5]; % Where the geodesics start [x;y]. You can set multiple seeds, as for the tips below.
input.exportValues=1;


output=MatlabHFM_RiemannExtra(input);

imagesc(squeeze(output.values(:,10,10,:)))
