% Copyright Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay

% This file demonstrates anisotropic fast marching with the "Half-Disk" model, in dimension 3+1.
% Three physical dimensions, and one extra isotropic dimension.

% If you have not compiled this model yet, please run the CompileMexHFM script and then
% compileModelsHFM(binary_Dir,JMM_CPPLibs_Dir,{'AsymmetricQuadratic3p1'})


% ------ 
% The metric applied to a vector v in R^4 takes the form 
% F(v)^2 =  sqrt(F_3(v_3)^2+ bundleCost^2*v_4^2),
% where F_3 is a three dimensional asymmetric quadratic metric, applied to the 
% first three components of f, and the bundleCost is a scalar multiplying the 
% fourth component of v.

% The three dimensional metric F_3 can be defined in three ways : 
% - Half disk model
% - Primal definition
% - Dual definition
% See the AsymmetricQuadratic2 demo for details.

clear input;

nx=20;
modelName = 'AsymmetricQuadratic3p1';
input.sndOrder=1;
input.dims = [nx;nx+1;nx+2;nx+3];
input.exportValues=1;

% Description of the three first dimensions
input.origin=[0;0;0]; % grid position parameters
input.gridScale=1/nx;

% Description of the fourth dimension.
input.bundleScale = 1./nx;

nGeo=8;

% ------ Array ordering -----

% We use in this example the ndgrid convention (as opposed to the meshgrid convention),
% or ordering the axes as x,y,z,t (as opposed to y,x,z,t).
input.arrayOrdering = 'ColumnMajor';

% ----- We choose a half disk model. ------

% The front propagates :
% - at speed |omega| in the direction of omega.
% - at speed r|omega| in the directions orthogonal to omega. (
% - at speed epsRev*|omega| in the direction -omega.
% - at speed 1/bundleCost in the fourth dimension.

% Parameters omega,r,s are provided by the user.
% Relaxation paramter epsRev as well.

input.halfDisk = zeros([3+1,nx,nx+1,nx+2,nx+3]); 

input.halfDisk(1,:,:,:,:)=1;      % omega_x
input.halfDisk(2,:,:,:,:)=0;      % omega_y
input.halfDisk(3,:,:,:,:)=0;      % omega_z
input.halfDisk(4,:,:,:,:)=1;      % r


input.bundleCost(:,:,:,:)=ones([nx,nx+1,nx+2,nx+3]);    % bundle inverse speed

input.seeds=[0.5;0.5;0.5;0.5]; % Starting point for the front propagation.


[x,y,z,t]=ndgrid(1:nx,1:(nx+1),1:(nx+2),1:(nx+3));
h=input.gridScale;
x=(x-0.5)*h;y=(y-0.5)*h;z=(z-0.5)*h;t=(t-0.5)*h;
% Center relative to the seed
x=x-0.5;y=y-0.5;z=z-0.5;t=t-0.5;
epsRev=0.2;

exactSolution = sqrt(max(x,0).^2+min(x,0).^2/epsRev.^2+y.^2+z.^2+t.^2);


output=MatlabHFM_AsymmetricQuadratic3p1(input);
imagesc(squeeze(output.values(:,10,10,:)))
pause;

diff=output.values-exactSolution;
imagesc(squeeze(diff(:,10,10,:)))

disp('largest error');disp(max(diff(:)));disp(min(diff(:)) );
