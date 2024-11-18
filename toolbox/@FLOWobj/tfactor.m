function [MS,CCL] = tfactor(FD,S,options)

%TFACTOR Transverse topographic symmetry (T-)factor
%
% Syntax
%
%     [MS,CL] = tfactor(FD,S,pn,pv,...)
%
% Description
%
%     tfactor calculates the transverse topographic symmetry factor also
%     called T-factor.  The T-factor quantifies the drainage basin
%     asymmetry and is calculated from the distance of the main river of a
%     drainage basin to the basin midline normalized by the distance of the
%     basin midline to the divide. The algorithm is explained in this blog
%     post https://topotoolbox.wordpress.com/2022/02/19/calculating-the-transverse-topographic-symmetry-t-factor/
%
%     Note that the function does not any correction along the DEM edges.
%     Drainage basins located on DEM margins are likely to feature wrong
%     t-factors. Check the function STREAMobj/removeedgeeffects.
%
%     Also note that the centerline may not be correct if a drainage
%     basin has a very irregular shape.
%
% Input arguments
%
%     FD     FLOWobj
%     S      STREAMobj of trunk rivers 
%     
%     Parameter name/value pairs
%     'K'    Smoothing factor (see STREAMobj/smooth) for smoothing the
%            planform of the centerline
%     'plot' {false} or true. If true, the function returns a map with
%            drainage basins colored with the tfactor, the centerlines and
%            the trunk rivers.
%     'waitbar' {true} or false
%     'geotable' {true} or false. If true, MS and CL will be geotables.
%            Otherwise, the function will return map structures. geotables
%            contain referencing information and thus can be plotted using
%            geoplot and exported as shapefiles including a .prj-file.
%
% Output arguments
%
%     MS     Mapping structure (or geotables) of drainage basins including 
%            the field 'tfactor'. MS contains the following fields:
%            .Geometry - Geometry type required for export to shapefile
%            .X, .Y    - Coordinates
%            .areapx   - Area in pixels
%            .Dd       - Distance to divides
%            .Da       - Distance to stream
%            .tfactor  - Da/Dd
%     CL     Structure array (or geotable) with centerlines. CL contains 
%            the following fields (or variables):
%            .Geometry - Geometry type required for export to shapefile
%            .X, .Y    - Coordinates
%            .DistToDivide - Distance to divide from centerline
%            .DistToStream - Distance to stream
%            .SmoothFactor - Smoothing factor (see STREAMobj/smooth) used
%            to smooth the centerline. 
%
% Example
%
%     DEM = readexample('taiwan');
%     DEM = inpaintnans(DEM);
%     S = trunk(STREAMobj(FD,'minarea',1000));
%     [MS,CL] = tfactor(FD,S,'plot',true);
%
%     % Both outputs can be exported to shapefiles
%     shapewrite(MS,'basins.shp')
%     shapewrite(CL,'centerlines.shp')
%
%
% See also: dbasymmetry, STREAMobj/removeedgeeffects, STREAMobj/smooth
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 15. November, 2024

arguments
    FD FLOWobj
    S  STREAMobj
    options.K  (1,1) {mustBeNumeric,mustBeNonnegative} = 0
    options.plot (1,1) = false
    options.waitbar (1,1) = false
    options.geotable (1,1) = true
end


if ~istrunk(S)
    error('TopoToolbox:wrongInput','S must be a network of trunk rivers.')
end

K = options.K;

%%
CS  = STREAMobj2cell(S);
xyO = cellfun(@(S) streampoi(S,'outlet','xy'),CS,'UniformOutput',false);
xyO = vertcat(xyO{:});
IXO = coord2ind(FD,xyO(:,1),xyO(:,2));

DB  = drainagebasins(FD,IXO);
nrbasins = double(max(DB));

% Preallocate 
% Structure array for midline
CCL  = struct('Geometry',repmat({'Line'},nrbasins,1),...
             'ID',num2cell(1:nrbasins)', ...
             'X',0,...
             'Y',0,...       
             'DistToDivide',0, ...
             'DistToStream',0, ...
             'SmoothFactor',0);
% Structure array for basin outlines
MS  = struct('Geometry',repmat({'Polygon'},nrbasins,1),...
             'tfactor',0, ...
             'areapx',0);
wb  = options.waitbar;
if wb
h   = waitbar(0,'Please wait'); 
end

for r = 1:nrbasins
    if wb
    waitbar(r/nrbasins,h)
    end
    
    % Identify basin
    IDB = DB == r;
    
    IDB = pad(crop(IDB,IDB),1);
    
    % Calculate first distance transform
    D   = GRIDobj(IDB);
    D.Z = bwdist(~IDB.Z,'e');
    
    % Calculate second distance transform
    D = 1./D;
    outletix = coord2ind(IDB,xyO(r,1),xyO(r,2));
    D.Z = graydist(D.Z,outletix,'q');
    D.Z(isinf(D.Z)) = nan;
    
    % Get centerline
    FDCL   = FLOWobj(D,'preprocess','carve');
    DDIST  = flowdistance(FDCL);
    
    [~,ix] = max(DDIST.Z(:));
    CL   = STREAMobj(FDCL,'chan',ix);

    % Smooth centerline
    if options.K ~= 0
        x = smooth(CL,CL.x,'K',K);
        y = smooth(CL,CL.y,'K',K);
    else
        x = CL.x;
        y = CL.y;
    end
    
    [~,~,x,y] = STREAMobj2XY(CL,x,y);
    
    % Remove points outside drainage basin which might occur due to
    % smoothing
    IX = coord2ind(IDB,x,y);
    inan = isnan(IX);
    IX(inan) = [];
    x(inan)  = [];
    y(inan)  = [];
    
    % Calculate distances
    [DIST,L]  = distance(IDB,IX);
    
    DBPERIM = bwperim(IDB.Z);
    dd = accumarray(L(DBPERIM),DIST.Z(DBPERIM),...
                    [numel(IDB.Z) 1],@max,...
                    nan(1,class(DIST.Z)));
    Dd = dd(IX);

    %% Deflection of the trunk stream from the centerline
    IXS = coord2ind(IDB,CS{r}.x,CS{r}.y);
    da = accumarray(L(IXS),DIST.Z(IXS),...
                    [numel(IDB.Z) 1],@max,...
                    nan(1,class(DIST.Z)));
    Da = da(IX);
    
    % Write centerline data
    CCL(r).X = x;
    CCL(r).Y = y;
    CCL(r).DistToDivide = Dd;
    CCL(r).DistToStream = Da;
    CCL(r).SmoothFactor = K;
    
    % Calculate T
    Da = mean(Da,'omitnan');
    Dd = mean(Dd,'omitnan');
    T  = Da./Dd;
    
    % Write data to structure arrays
    IMS = GRIDobj2polygon(IDB);
    
    MS(r).X = IMS.X;
    MS(r).Y = IMS.Y;
    MS(r).tfactor = double(T);
    MS(r).Da = Da;
    MS(r).Dd = Dd;
    MS(r).areapx = nnz(IDB.Z);
    
end

if wb
    close(h)
end

if options.geotable
    MS  = mapstruct2geotable(MS,"coordinateSystemType","planar",...
        "CoordinateReferenceSystem",parseCRS(S));
    CCL =  mapstruct2geotable(CCL,"coordinateSystemType","planar",...
        "CoordinateReferenceSystem",parseCRS(S));
end

if options.plot
    figure
    
    colorRange = makesymbolspec('Polygon',...
                           {'tfactor',[0 1],'FaceColor',summer(10)});
    mapshow(MS,'SymbolSpec',colorRange)
    hold on
    plot(S,'color','b')
    mapshow(CCL,'color','k')
    
    colormap(summer(10))
    clim([0 1])
    h = colorbar;
    h.Label.String = 'T-factor';
end
