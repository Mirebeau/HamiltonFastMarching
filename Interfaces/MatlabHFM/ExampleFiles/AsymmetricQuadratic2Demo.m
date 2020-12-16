% Copyright Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay

% This file demonstrate anisotropic fast marching with a metric of 
% AsymmetricQuadratic type, in dimension two.
% Higher dimensional counterparts are also available, in dimension <= 5.

% If you have not compiled this model yet, please run the CompileMexHFM script and then
% compileModelsHFM(binary_Dir,JMM_CPPLibs_Dir,{'AsymmetricQuadratic2'})

% A metric of AsymmetricQuadratic type takes the form
% F(v)^2 = <v,M v> + <omega,v>_+^2
% The symmetric matrix M and vector omega are parameters that may depend on
% the current point. The matrix M must be positive definite.

% The dual metric, defined in the co-tangent space, has a similar form.
% F^*(v)^2 = <v,D v> + <eta,v>_+^2.

% This type of metric can be used to approximate 'HalfDisk models', whose unit ball 
% is a half disk, or half ellipse, up to a small relaxation parameter.

% There are three ways to specify the metric:
% - directly define M and omega
% - define the dual parameters D and eta
% - define a half disk model.

clear input;
nx=81;ny=83; %Taking different dimensions to make sure axes are not inverted...
nGeo=8;

modelName = 'AsymmetricQuadratic2';
metricType = 'Dual'; % Primal % Dual

% By default, we do use Matlab's' convention of ordering the axes as y, x, z, 
% (and relying on ColumnMajor array data indexation).
% use input.arrayOrdering = 'ColumnMajor' for standard axes ordering.

input.origin=[0;0]; % grid position parameters
input.gridScale=1/nx;
input.dims = [nx;ny];

input.seeds=[0.5;0.5]; % Seed for the front propagation.

% Create the corresponding coordinate system
[x,y] = meshgrid(1:nx,1:ny);
x=(x-0.5)*input.gridScale; y=(y-0.5)*input.gridScale;
X = x-0.5; Y=y-0.5;
% ------------- Half Disk model ------------
if strcmp(metricType,'HalfDisk')
    % The front goes in the direction of a given vector 
    % omega, at speed |omega| (the norm of omega),
    % and in the sideways direction at speed r |omega|.
    % Main inputs are omega, r, and epsRev (relaxation parameter).
    % Motion in the direction of -omega is in principle forbidden 
    % (in practice, it still happens at speed epsRev*r*|omega|).
   
    % The hamiltonian is F^*(v) = <omega,v>_+^2 + ...
    % r^2 |omega|^2 |P v|^2 + ...
    % r^2 epsRev^2 <omega,v>_+^2,
    % where omega=(omega_x,omega_y), and P is the orthogonal projector 
    % onto the plane orthogonal to omega.
    
    input.halfDisk = zeros([3,ny,nx]);
    input.halfDisk(1,:,:)=1;      % omega_x
    input.halfDisk(2,:,:)=0.5;    % omega_y
    input.halfDisk(3,:,:)=0.5;    % r

    % Relaxation parameters are:
    % - a scalar eps, for how accurately forward speed is implemented.
    % - a scalar epsRev, for how stongly backwards motion is penalized.
    % As usual, small parameters yield large stencils. We leave default values here.
    
    omega = [1.,0.5]';
    eta = omega/(omega'*omega);
    epsRev=0.2;
    r=0.5;
    exactSolution = sqrt( max(0,eta(1)*X+eta(2)*Y).^2+ ...
        min(0,eta(1)*X+eta(2)*Y).^2/(r*epsRev).^2+ ...
        (eta(2)*X-eta(1)*Y).^2/r^2);
    
% --------- Primal definition ---------
elseif strcmp(metricType,'Primal')
    input.metric = zeros([5,ny,nx]);
    input.metric(1,:,:) = 1; % Mxx
    input.metric(2,:,:) = 0; % Mxy
    input.metric(3,:,:) = 1; % Myy
    input.metric(4,:,:) = 1; % omega_x
    input.metric(5,:,:) = 1; % omega_y
    
    exactSolution = sqrt(X.^2+Y.^2+max(X+Y,0).^2);
    
% ---------- Dual definition ---------
elseif strcmp(metricType,'Dual')
    input.dualMetric = zeros([5,ny,nx]);
    input.dualMetric(1,:,:) = 1; % Dxx
    input.dualMetric(2,:,:) = 0; % Dxy
    input.dualMetric(3,:,:) = 1; % Dyy
    input.dualMetric(4,:,:) = 1; % eta_x
    input.dualMetric(5,:,:) = 1; % eta_y
    
    exactSolution = sqrt(X.^2+Y.^2-max(X+Y,0).^2/3);
    
else 
    assert(false,'unrecognized metric type')
end

%        input.stopWhenTipsReached=1;

% We set multiple tips.
epsGeo = 1/nGeo;
arrGeo = (epsGeo/2):epsGeo:(1-epsGeo/2);
[x,y]=meshgrid(arrGeo,arrGeo);
x=reshape(x,[numel(x),1]); y=reshape(y,[numel(y),1]);
input.tips = [x,y]'; % where the geodesics end. This will select the angle t for which the geodesic is the shortest.
input.exportValues=1; % distance table, of size [n,n] (minimum of the previous one over directions)
input.order=2;
output=eval(['MatlabHFM_' modelName '(input)']);

clf;
imagesc(output.values)
geodesics = mat2cell(output.geodesicPoints,2,output.geodesicLengths);
for i=1:size(geodesics,2) %geodesics joining the "tips"
    rescaledGeodesic=RescaledCoords(geodesics{i},input.origin,[input.gridScale;input.gridScale]);
    line(rescaledGeodesic(1,:),rescaledGeodesic(2,:));
end

%imagesc(output.values-exactSolution) %Control
