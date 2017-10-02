if true %Some geodesics, around the origin, without obstacles.
    clear input;
    n=101;
    nTheta=60;
    nGeo=8;
    
    % Defining the domain
    input.dims = [n;n;nTheta];
    input.origin=[0;0];  % Physical origin
    input.gridScale=1/n; % Physical gridScale

    input.model = 'Elastica2<5>';%Alternatively : 'ReedsShepp2','ReedsSheppForward2', 'Elastica2<5>', 'Dubins2'
    input.speed = 1; %Use an array for a position, angle, or (position and angle) dependent speed. 
    % Alternatively input.speed = ones(n,n), or ones(nTheta), or ones(n,n,nTheta)
    input.xi = 0.1; %Model parameter, typical radius of curvature.
    input.eps=0.1;
    
    input.seeds=[0.5;0.5;0]; % Where the geodesics start [x;y;theta]. You can set multiple seeds, as for the tips below.
    % Defining the geodesic endpoints
    epsGeo = 1/nGeo;
    arrGeo = (epsGeo/2):epsGeo:(1-epsGeo/2);
    [x,y]=meshgrid(arrGeo,arrGeo);
    x=reshape(x,[numel(x),1]); y=reshape(y,[numel(y),1]);
    input.tips = [x,y,(pi/3)*ones(size(x))]'; % where the geodesics end. [x;y;theta]
    input.exportValues=1; % distance table, of size [n,n,numberOfDirections]
    
    input.pointToIndex=input.seeds;
    output=MatlabHFM_Curvature2(input);
    
    clf;
    imagesc(min(output.values,[],3))
    geodesics = mat2cell(output.geodesicPoints,3,output.geodesicLengths);    
    for i=1:size(geodesics,2) %geodesics joining the "tips"
        rescaledGeodesic=RescaledCoords(geodesics{i}(1:2,:),input.origin,[input.gridScale;input.gridScale]);
        line(rescaledGeodesic(1,:),rescaledGeodesic(2,:));
    end;
    pause;
end

if false % Within centre Pompidou
    clear input;
    input.model = 'ReedsSheppForward2'; % Alternatively 'ReedSheppForward2', %'Elastica2<5>', 'Dubins2'; 
    input.xi = 0.7; %Model parameter, typical radius of curvature.
    input.eps=0.1;
    input.speed=1;
    input.projective=1; % Applies to the ReedsShepp model only
    
    im = imread('centre_pompidou_800x546.png');
    im = ( (im(:,:,1)==255) | (im(:,:,3)==255)) & (im(:,:,2)==0);    
    input.walls = 1.-double(im);
        
    gridScale=1./90;
    input.origin=[0;0];
    input.gridScale=gridScale;
    input.seeds_Unoriented=[80,80;170,290]*gridScale;
    
    input.exportValues=1;
    input.dims = [fliplr(size(im(:,:,1)))';60/(1+input.projective)];
    
    input.tips_Unoriented=...
        [369.4, 252.2, 285., 418.6, 479.8, 687.2, 745.8, 740.4, 593.8, 558.6,...
        599.2, 497.2, 495.8, 427.2, 339., 264.6, 242.4, 354.6, 191.6, ...
        178.8, 105.8, 124., 127., 419.2; ...
        482.5, 354.5, 478., 488., 487.5, ...
        478., 502.5, 300., 225.5, 378., 475.5, 81., 127.5, 128., 111., 108.,...
        176.5, 290.5, 110., 252.5, 428.5, 494., 353., 421.] *gridScale;    
    input.geodesicVolumeBound=12;
    output=MatlabHFM_Curvature2(input);
    
    dist=output.values;
    dist(dist==Inf)=0;
    
    clf;
    imagesc(min(dist,[],3))
    geodesics = mat2cell(output.geodesicPoints_Unoriented,3,output.geodesicLengths_Unoriented);    
    for i=1:size(geodesics,2)
        rescaledGeodesic=RescaledCoords(geodesics{i}(1:2,:),input.origin,[input.gridScale;input.gridScale]);
        line(rescaledGeodesic(1,:),rescaledGeodesic(2,:));
    end;

end

if false % Straight geodesics, for comparison
    clear input;
    input.model = 'IsotropicBox2<Boundary::Closed>'; % Alternatively 'ReedSheppForward2', %'Elastica2<5>', 'Dubins2'; 
    input.speed=1;
    
    im = imread('centre_pompidou_800x546.png');
    im = ( (im(:,:,1)==255) | (im(:,:,3)==255)) & (im(:,:,2)==0);    
    input.walls = 1.-double(im);
        
    gridScale=1./90;
    input.origin=[0;0];
    input.gridScale=gridScale;
    input.seeds=[80,80;170,290]*gridScale;
    
    input.exportValues=1;
    input.dims = fliplr(size(im(:,:,1)))';
    
    input.tips=...
        [369.4, 252.2, 285., 418.6, 479.8, 687.2, 745.8, 740.4, 593.8, 558.6,...
        599.2, 497.2, 495.8, 427.2, 339., 264.6, 242.4, 354.6, 191.6, ...
        178.8, 105.8, 124., 127., 419.2; ...
        482.5, 354.5, 478., 488., 487.5, ...
        478., 502.5, 300., 225.5, 378., 475.5, 81., 127.5, 128., 111., 108.,...
        176.5, 290.5, 110., 252.5, 428.5, 494., 353., 421.] *gridScale;    
    
    output=MatlabHFM_Isotropic(input);
    
    dist=output.values;
    dist(dist==Inf)=0;
    
    clf;
    imagesc(dist)
    geodesics = mat2cell(output.geodesicPoints,2,output.geodesicLengths);    
    for i=1:size(geodesics,2)
        rescaledGeodesic=RescaledCoords(geodesics{i}(1:2,:),input.origin,[input.gridScale;input.gridScale]);
        line(rescaledGeodesic(1,:),rescaledGeodesic(2,:));
    end;

end
