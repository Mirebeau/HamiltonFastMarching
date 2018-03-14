% If you have not compiled this model yet, please run the CompileMexHFM script and then
% compileModelsHFM(binary_Dir,JMM_CPPLibs_Dir,{'RiemannLifted2_Periodic'})

if true 
    % This example reproduces the the Reeds Shepp model.
    clear('input')
    n=101;
    nTheta=60;
    input.dims = [n;n;nTheta];
    input.gridScale = 1./n;
    input.bundleScale = 2.*pi/nTheta; %This is actually the default value.
    input.origin=[0;0]; %Physical origin. This is actually the default value.
    
    
    % We choose a (2+1)D model, where the first two dimensions are 
    % equipped with a Riemannian metric, and the last dimension is periodic.
    modelName = 'RiemannLifted2_Periodic';
    % If you want a non-periodic last dimension, then select 'RiemannLifted2<Boundary::Closed>'
    
     
    %The dual metric is the inverse of the metric.
    % In this model, it is block diagonal, with structure
    % mxx, mxy, myy, mzz    
    input.dualMetric = zeros([4,n,n,nTheta]);
    
    % The physical part of the metric is strongly anisotropic, 
    % and aligned with the angle theta
    eps = 0.1;
    for i=1:nTheta
        theta = (2.*pi*(i-1))/nTheta;
        input.dualMetric(1,:,:,i) = cos(theta)^2+eps^2*sin(theta)^2;
        input.dualMetric(2,:,:,i) = cos(theta)*sin(theta)*(1-eps^2);
        input.dualMetric(3,:,:,i) = sin(theta)^2+eps^2*cos(theta)^2;
    end
    
    % The bundle part of the metric is constant
    xi=0.3; %Homogeneous to a radius of curvature
    input.dualMetric(4,:,:,:) = 1./xi^2;

    input.seeds = [0.5;0.5;0.];
    
    nGeo=8;
    epsGeo = 1/nGeo;
    arrGeo = (epsGeo/2):epsGeo:(1-epsGeo/2);
    [x,y]=meshgrid(arrGeo,arrGeo);
    x=reshape(x,[numel(x),1]); y=reshape(y,[numel(y),1]);
    input.tips = [x,y,(pi/3)*ones(size(x))]'; % where the geodesics end. [x;y;theta]

%    input.tips = [0.2,0.7,0.8;0.6,0.2,0.5;pi/3,pi/2,-pi];
    input.exportValues=1; 
%    input.pointToIndex = input.seeds;
    output = eval(['MatlabHFM_' modelName '(input)']);
    
    clf;
    imagesc(min(output.values,[],3))
    geodesics = mat2cell(output.geodesicPoints,3,output.geodesicLengths);    
    for i=1:size(geodesics,2) %geodesics joining the "tips"
        rescaledGeodesic=RescaledCoords(geodesics{i}(1:2,:),input.origin,[input.gridScale;input.gridScale]);
        line(rescaledGeodesic(1,:),rescaledGeodesic(2,:));
    end
end