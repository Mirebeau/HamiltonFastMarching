% Demo of PrescribedCurvature2 executable. 
% It must be compiled first by executing CompileHFM(binary_Dir,JMM_CPPLibs_Dir,'PrescribedCurvature2'), from the CompileHFM.m script.
% The cost takes the form C(xi(x''-kappa))/speed, where speed, xi and kappa are either constant or state dependent.
% The function C is model dependent, see the Curvature2Demo script for details.

if true %Some geodesics, around the origin, without obstacles.
    clear input;
    n=101;
    nTheta=60;
    nGeo=8;
    
    % Defining the domain
    input.dims = [n;n;nTheta];
    input.origin=[0;0];  % Physical origin
    input.gridScale=1/n; % Physical gridScale
    arrPos = (input.gridScale/2):input.gridScale:(1-input.gridScale/2);
    [x,y,theta]=meshgrid(arrPos,arrPos,(0:(nTheta-1))*(2*pi/nTheta));
    
    modelName = 'DubinsExt2'; %Alternatively : 'ReedsSheppExt2','ReedsSheppForwardExt2', 'ElasticaExt2', 'DubinsExt2'
    input.eps=0.1;%Relaxation parameter

    % The three folloinwg fields may be uniform over the domain, or state dependent. E.g. =1, or =ones(n,n,nTheta).
    
 %   input.speed = 1; input.xi = 0.15; input.kappa=1./(2*input.xi); % A car that turns better left than right.
    input.speed = 1; input.xi = 0.15./(1+(y<0.5)); input.kappa=0; % A car that turns better on the upper part of the domain.
    
        % Not currently available :
    % - Fields depending only on position, e.g. ones(n,n), or on orientation, e.g. ones(nTheta). (Could be if requested)
    % - Time dependency and automatic differentiation (could be for speed if requested)
    
    input.seeds=[0.5;0.5;0]; % Where the geodesics start [x;y;theta]. You can set multiple seeds, as for the tips below.
    % Defining the geodesic endpoints
    epsGeo = 1/nGeo;
    arrGeo = (epsGeo/2):epsGeo:(1-epsGeo/2);
    [x,y]=meshgrid(arrGeo,arrGeo);
    x=reshape(x,[numel(x),1]); y=reshape(y,[numel(y),1]);
    input.tips = [x,y,(pi/3)*ones(size(x))]'; % where the geodesics end. [x;y;theta]
    input.exportValues=1; % distance table, of size [n,n,numberOfDirections]
    
    input.pointToIndex=input.seeds;
    output=eval(['MatlabHFM_' modelName '(input)']);
    
    clf;
    imagesc(min(output.values,[],3))
    geodesics = mat2cell(output.geodesicPoints,3,output.geodesicLengths);    
    for i=1:size(geodesics,2) %geodesics joining the "tips"
        rescaledGeodesic=RescaledCoords(geodesics{i}(1:2,:),input.origin,[input.gridScale;input.gridScale]);
        line(rescaledGeodesic(1,:),rescaledGeodesic(2,:));
    end
end