function DEMr = project(SOURCE,TARGET,varargin)

%PROJECT Warps a GRIDobj to a different coordinate system
%
% Syntax
%
%     OUT = project(SOURCE,TARGET)
%     OUT = project(SOURCE,TARGET,pn,pv,...)
%
% Description
%
%     project reprojects the GRIDobj 'SOURCE' to have the same projection
%     as the GRIDobj 'TARGET'. Both GRIDobj's need to have a valid
%     projection in their 'georef' properties. If SOURCE does not, then it
%     is assumed to have a geographic coordinate system (horizontal WGS84
%     datum). The function supports to project from geographic to projected
%     coordinates, and from projected to geographic coordinates. It also
%     supports transformation from a projected to another projected
%     coordinate (via WGS84 geographic coordinates). The spatial resolution
%     of the output (OUT) will be the same as of TARGET, unless set
%     differently by the optional variable res.
%
% Input arguments
%
%     SOURCE         instance of GRIDobj that shall be transformed
%     TARGET         instance of GRIDobj with the target projection. 
%                    Alternatively, target can be a raster reference object
%                    (MapCellReference, GeoCellReference, ...) or it can be
%                    a EPSG number which is accepted by the function
%                    projcrs or geocrs. Note that there are also projection
%                    IDs from the authority ESRI. In this case, use
%                    projcrs(number,'authority','ESRI').
%
%     Parameter name/value pairs
%     
%     'res'          scalar
%                    spatial resolution. Default is TARGET.cellsize. If
%                    target is an mstruct, res must be supplied. 
%     'method'       string
%                    'bilinear' (default),'bicubic' or 'nearest'
%     'fillvalue'    scalar
%                    nan (for single and double grids). Otherwise 0.
%     'align'        logical scalar
%                    true (default) or false. If true, the OUT grid will be
%                    spatially aligned with TARGET. This means, they have
%                    the same extent and spatial alignment of cells. If
%                    'align', true then setting the resolution 'res' will
%                    have no effect.
%                    
% Output arguments
%
%     OUT            instance of GRIDobj
%
% Example
%
%    DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%    DEMgcs = project(DEM,4326,res = 0.001);
%    imageschs(DEMgcs)
%
% See also: GRIDobj, GRIDobj/reproject2utm, egm96heights, imtransform
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de) and Wolfgang
%         Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 31. May, 2024

% Check input arguments
narginchk(2,inf)

[cs,prjtarget,targetisGRIDobj,targetisproj] = checktarget(TARGET);

p = inputParser;
p.FunctionName = 'GRIDobj/project';
addParameter(p,'res', cs,@(x) isscalar(x));
addParameter(p,'align',true,@(x) isscalar(x));
addParameter(p,'method','linear',@(x) ischar(x));
addParameter(p,'fillvalue',nan,@(x) isscalar(x));
addParameter(p,'imtransform',true)
parse(p,varargin{:});

[prjsource,sourceisproj] = checksource(SOURCE);

% Determine transformation type
pr2pr  = sourceisproj && targetisproj;
pr2gcs = sourceisproj && ~targetisproj;
gcs2pr = ~sourceisproj && targetisproj;
gcs2gcs = ~sourceisproj && ~targetisproj;
if gcs2gcs
    error("Geographic-to-geographic coordinate projection is not supported.")
end
% Make sure that only one option is chosen
if nnz([pr2pr pr2gcs gcs2pr]) ~= 1
    error('Something went wrong.')
end

% Check, if cellsize was supplied
if ~targetisGRIDobj
    cs = p.Results.res;
    if isempty(cs)
        error('Spatial resolution of the target GRIDobj must be set. See option "res".');
    end
end
        

% First, we create the reference geometry of the target grid
if ~p.Results.align || ~targetisGRIDobj 
    % If we only know the target projection, then we need to find the
    % boundaries of the source and project them to the target projection
    [x0,y0]   = getcoordinates(SOURCE);
    x0        = x0(:);
    y0        = y0(:);

    % Calculate bounds of the reprojected DEM
    x         = [x0; x0; repmat(x0(1),numel(y0),1); repmat(x0(end),numel(y0),1)];
    y         = [repmat(y0(1),numel(x0),1); repmat(y0(end),numel(x0),1); y0; y0]; 
    
    if pr2pr
        [lat,lon] = projinv(prjsource,x,y);
        [x,y]     = projfwd(prjtarget,lat,lon);
    elseif gcs2pr
        [x,y]     = projfwd(prjtarget,y,x);
    elseif pr2gcs
        [y,x]     = projinv(prjsource,x,y);
    else
        error("Something's wrong here...")
    end
    
    % We should actually use imref2d directly. Since we need to set the
    % cell extents precisely, we use maprefcells as an intermediate step
    
    % We assume that the target grid is cell referenced
    xextent = [min(x),max(x)];
    yextent = [min(y),max(y)];
    if ~pr2gcs
        R = maprefcells(xextent,yextent,cs,cs);
        Rtarget = imref2d(R.RasterSize,R.XWorldLimits,R.YWorldLimits);
    else
        R = georefcells(yextent,xextent,cs,cs);
        Rtarget = imref2d(R.RasterSize,R.LongitudeLimits,R.LatitudeLimits);
    end

else
    % If we have a target grid, then it is easy
    Rtarget = GRIDobj2imref2d(TARGET,pr2pr || gcs2pr); 
end

% Limits of source grid
[Rsource,Zsource] = GRIDobj2imref2d(SOURCE,pr2pr || pr2gcs); 

% get fillvalue
fillval = p.Results.fillvalue;
if isinteger(SOURCE.Z)
    fillval = zeros(1,class(SOURCE.Z));
elseif islogical(SOURCE.Z)
    fillval = false;
else
    fillval = cast(fillval,class(SOURCE.Z));
end
fillval = double(fillval);
    
% check method
if p.Results.imtransform
    meth = p.Results.method;
    switch meth
        case {'linear','cubic'}
            meth = ['bi' char(meth)];
    end
    meth = validatestring(meth,{'bilinear','nearest','bicubic'});
else
    % if imwarp is use
    meth = validatestring(p.Results.method,{'linear','nearest','cubic'});
end

% Prepare tform for the image transform
if pr2pr
    if ~p.Results.imtransform
        T = geometricTransform2d(@FWDTRANSpr2pr, @INVTRANSpr2pr);
    else
        T = maketform('custom', 2, 2, @FWDTRANSpr2pr, @INVTRANSpr2pr, []);
    end
elseif gcs2pr
    if ~p.Results.imtransform
        T = geometricTransform2d(@FWDTRANSgcs2pr, @INVTRANSgcs2pr);
    else
        T = maketform('custom', 2, 2, @FWDTRANSgcs2pr, @INVTRANSgcs2pr, []);
    end
else % pr2gcs
    if ~p.Results.imtransform
        T = geometricTransform2d(@FWDTRANSpr2gcs, @INVTRANSpr2gcs);
    else
        T = maketform('custom', 2, 2, @FWDTRANSpr2gcs, @INVTRANSpr2gcs, []);
    end
end

% Calculate image transform 
if ~p.Results.imtransform
    Znew = imwarp(Zsource,Rsource,T,meth,...
        'OutputView',Rtarget,'FillValues',fillval);
else

    Znew = imtransform(Zsource,T,meth,...
        'Xdata',Rtarget.XWorldLimits,'Ydata',Rtarget.YWorldLimits,...
        'Udata',Rsource.XWorldLimits,'Vdata',Rsource.YWorldLimits,...
        'XYScale',[Rtarget.PixelExtentInWorldX Rtarget.PixelExtentInWorldY],...
        'Fillvalues',fillval,'size',Rtarget.ImageSize); % ,xdata,ydata
end

if p.Results.align && targetisGRIDobj
    % If aligned, this is easy. The transformed GRIDobj will be perfectly
    % aligned with TARGET
    DEMr = GRIDobj(TARGET,Znew);
    % DEMr.Z = flipud(DEMr.Z);
else
    % We have calculated the imtransform with 'ColumnsStartFrom' south. 
    % GRIDobjs use 'ColumnsStartFrom' north
    DEMr = GRIDobj([]);
    DEMr.Z = Znew;
    DEMr.cellsize   = cs;
    DEMr.size = size(Znew);
    DEMr.georef = R;
    DEMr.wf  = worldFileMatrix(R);
    if pr2pr || gcs2pr
        DEMr.georef.ProjectedCRS = prjtarget;
    else
        DEMr.georef.GeographicCRS = prjtarget;
    end
end

DEMr.name = [SOURCE.name ' (projected)'];

% test output for equal sizes
if ~isequal(DEMr.size,DEMr.georef.RasterSize)
    error("Inconsistent raster size of georef and grid")
end
% function ends here

%% -----------------------------------------------------------------------
% Transformation functions for imtransform (projected --> projected)
    function x = FWDTRANSpr2pr(u,~)
        [tlat,tlon] = projinv(prjsource,u(:,1),u(:,2));
        [tx,ty] = projfwd(prjtarget,tlat,tlon);
        x = [tx ty];        
    end

 
    function u = INVTRANSpr2pr(x,~)
        [tlat,tlon] = projinv(prjtarget,x(:,1),x(:,2));
        [tx,ty] = projfwd(prjsource,tlat,tlon);
        u = [tx ty];
    end

% Transformation functions for imtransform (gcs --> projected)
    function x = FWDTRANSgcs2pr(u,~)
        [tx,ty] = projfwd(prjtarget,u(:,2),u(:,1));
        x = [tx ty];        
    end

 
    function u = INVTRANSgcs2pr(x,~)
        [tlat,tlon] = projinv(prjtarget,x(:,1),x(:,2));
        u = [tlon tlat];
    end

% Transformation functions for imtransform (projected --> gcs)
    function x = FWDTRANSpr2gcs(u,~)
        [tlat,tlon] = projinv(prjsource,u(:,1),u(:,2));
        x = [tlon tlat];        
    end

    function u = INVTRANSpr2gcs(x,~)
        [tx,ty] = projfwd(prjsource,x(:,2),x(:,1));
        u = [tx ty];
    end

end


%% -----------------------------------------------------------------------
function [cs,prj,targetisGRIDobj,targetisproj] = checktarget(TARGET)
% Check source
if isa(TARGET,'GRIDobj')
    % target is a GRIDobj. If yes, there's a cellsize.
    cs = TARGET.cellsize;
    targetisGRIDobj = true;

    % The target grid should have a projcrs object. If not, then we assume
    % that the target grid is in geographic coordinates and we will do an
    % inverse projection.
    try
        prj = TARGET.georef.ProjectedCRS;
        targetisproj    = true;
        
    catch
        try
            prj = TARGET.georef.GeographicCRS;
            targetisproj    = false;
        catch
            prj = [];
            targetisproj    = false;
        end
        
    end
elseif isa(TARGET,"geocrs")
    % target is a geocrs object
    cs = [];
    targetisGRIDobj = false;
    targetisproj    = false;
    prj = TARGET;

elseif isa(TARGET,"projcrs")
    % target is a projcrs object
    cs = [];
    targetisGRIDobj = false;
    targetisproj    = true;
    prj = TARGET;

else
    % target is some epsg-number or some string that either projcrs or
    % geocrs can handle.
    cs = [];
    targetisGRIDobj = false;
    try
        prj = projcrs(TARGET);
        targetisproj = true;

    catch
        prj = geocrs(TARGET);
        targetisproj = false;
    end
end
end

%% -----------------------------------------------------------------------
function [prj,sourceisproj] = checksource(SOURCE)

% Get georef of SOURCE
if isempty(SOURCE.georef)
    sourceisproj = false;
    prj = [];
else
    if isProjected(SOURCE)
            prj = SOURCE.georef.ProjectedCRS;
            sourceisproj = true;
    else        
            prj = SOURCE.georef.GeographicCRS;
            sourceisproj = false;
    end
end
end

%% ----------------------------------------------------------------------
% functions to create imref2d objects
function [R,Z] = GRIDobj2imref2d(DEM,isproj)
RM = DEM.georef;
if isproj
    R = imref2d(RM.RasterSize,...
        RM.XWorldLimits,RM.YWorldLimits);
    columnsStartNorth = strcmp(RM.ColumnsStartFrom,'north');
else
    if ~isempty(RM)
        R = imref2d(RM.RasterSize,...
            RM.LongitudeLimits,RM.LatitudeLimits);
        columnsStartNorth = strcmp(RM.ColumnsStartFrom,'north');
    else
        [x,y] = getcoordinates(DEM);
        if y(1) > y(end)
            columnsStartNorth = true;
        end
        R = imref2d(DEM.size,[x(1) x(end)],[min(y) max(y)]);
    end
end
if nargout == 2
Z = DEM.Z;
if columnsStartNorth 
    Z = flipud(Z);
end
end
end


