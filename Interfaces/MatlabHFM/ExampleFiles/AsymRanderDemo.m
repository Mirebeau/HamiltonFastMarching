% Copyright Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay

% This file demonstrate anisotropic fast marching with 'asymmetric Rander metrics',
% a custom class of Finslerian metrics, available in two dimensions, and reading
% F(x) = sqrt(<x,Mx> + <x,u>_+^2 + <x,v>_+^2 ) + <x,w>
% where a_+ denotes max(a,0). The matrix m must be positive definite, the vectors u and v
% are arbitrary, the vector w must be sufficiently small so that F(x)>0 for all x!=0.

% The metric coordinates are specified as 
% m_xx, m_xy, m_yy, u_x, u_y, v_x, v_y, w_x, w_y 

% ADVERTISEMENT : A more extended demo is available in my python notebooks.
% It could easily be translated to Matlab. 
% github.com/Mirebeau/AdaptiveGridDiscretizations
% Fast Marching methods -> Asymmetric quadratic metrics

% If you have not compiled this model yet, please run the CompileMexHFM script and then
% compileModelsHFM(binary_Dir,{'AsymRander2'})

% We do use Matlab's convention of ordering the axes as y, x, z.
% Use the input.arrayOrdering field to change this behavior
clear hfmIn

% Generate a coordinate system
% In contrast with the other example files, we use here ndgrid instead of
% meshgrid. Thus coordinate axes are xy instead of yx.
h = 0.1;
[x,y]=ndgrid(-1:0.1:1,-1:0.1:0.5); % use xy coordinate axes instead of meshgrid's yx

% Specify it to the hfm library
hfmIn.gridScale=h;
hfmIn.origin = [x(1,1);y(1,1)] - h/2;
[nx,ny] = size(x);
hfmIn.dims = [nx;ny];
hfmIn.arrayOrdering = 'ColumnMajor'; % Adequate for ndgrid. (Remove for meshgrid.)

% Other input data
hfmIn.seeds = [0;0];
hfmIn.exportValues = 1;
%hfmIn.factoringRadius = 200; % Apply source factorization over 200 pixels radius

% Specify the metric: m_xx, m_xy, m_yy, u_x, u_y, v_x, v_y, w_x, w_y
hfmIn.metric = zeros([9,nx,ny]);

% First example : m=Id, u and v are axis aligned, w is zero.
% m=Id
hfmIn.metric(1,:,:)=1;
hfmIn.metric(2,:,:)=0;
hfmIn.metric(3,:,:)=1;
% u = [3;0]
hfmIn.metric(4,:,:)=3;
hfmIn.metric(5,:,:)=0;
% v = [0;3]
hfmIn.metric(6,:,:)=0;
hfmIn.metric(7,:,:)=3;
% w = [0;0]
hfmIn.metric(8,:,:)=0;
hfmIn.metric(9,:,:)=0;

hfmOut = MatlabHFM_AsymRander2(hfmIn);

% Since the metric is constant in space, the exact solution is known

numSolution = hfmOut.values;
exactSolution = sqrt(x.^2+y.^2+max(0,3*x).^2+max(0,3*y).^2);

max(abs(exactSolution - numSolution),[],'all') % max error

% Plot the values 
% (Todo : take into account that we use ndgrid, rather than meshgrid)
imagesc(numSolution) 
