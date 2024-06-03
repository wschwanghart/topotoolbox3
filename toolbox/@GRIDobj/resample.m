function DEMr = resample(DEM,target,options)

%RESAMPLE change spatial resolution of a GRIDobj
%
% Syntax
%
%     DEMr = resample(DEM,cellsize)
%     DEMr = resample(DEM,GRID)
%     DEMr = resample(...,method)
%     DEMr = resample(DEM,GRID,method,swapzone)
%
% Description
%
%     resample changes the cellsize of a grid. The function uses the MATLAB
%     function imtransform. If an instance of GRIDobj is supplied as
%     second argument, resample interpolates values in DEM to match the
%     spatial reference of GRID.
%
% Input arguments
%
%     DEM       grid object (GRIDobj)
%     cellsize  cellsize of resampled grid
%     GRID      other grid object
%     method    'bicubic', 'bilinear', or 'nearest' 
%     swapzone  true or false. If true and if DEM and GRID have different
%               projected coordinate systems, the function will attempt to
%               reproject and resample the DEM in one step. Note that this
%               requires the mapping toolbox. In case, the DEM is in a
%               geographic coordinate system, please use the function
%               reproject2utm(DEM,GRID).
%
% Output arguments
%
%     DEMr    grid object (GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEMr = resample(DEM,100);
%     imagesc(DEMr)
%
%
% See also: griddedInterpolant, imtransform
%        
% Author:  Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 1. June, 2024 

% check also mapresize, georesize

arguments
    DEM   GRIDobj
    target  {mustBeGRIDobjOrPositiveScalar}
    options.method   = "linear"
    options.swapzone = false
    options.crop     = true
end

method = validatestring(options.method,{'cubic', 'linear', 'nearest' },...
    'GRIDobj/resample','method',3);

% Fillvalues and 
if isUnderlyingType(DEM,'logical') || isUnderlyingInteger(DEM)
    fillval = 0;
    method = 'nearest';
else
    fillval = nan;
end

% get coordinate vectors
[Rsource,Zsource] = GRIDobj2imref2d(DEM,true);

% tform
% T = maketform('affine',);
T = affinetform2d([1 0 0; 0 1 0; 0 0 1]);

if isa(target,'GRIDobj')
    % the target is another GRIDobj
    DEMr    = target;
    % xlimitsn = DEMr.georef.XWorldLimits;
    % ylimitsn = DEMr.georef.YWorldLimits;
    % csxnew    = DEMr.georef.CellExtentInWorldX;
    % csynew    = DEMr.georef.CellExtentInWorldY;
    sizenew = DEMr.georef.RasterSize;
    Rtarget = GRIDobj2imref2d(DEMr,true);

    DEMr   = imwarp(Zsource,Rsource,T,method,...
        "FillValues",fillval,"OutputView",Rtarget);
    
    % DEMr.Z = imtransform(Zsource,T,method,...
    %     'Udata',[u(1) u(end)],'Vdata',[v(1) v(end)],...
    %     'Xdata',xlimitsn,'Ydata',ylimitsn,...
    %     'XYscale',[csxnew csynew],...
    %     'FillValues',fillval);
    if isequal(size(DEMr.Z),sizenew)
        error('Size after transformation inconsistent.')
    end
    DEMr.name = [DEM.name ' (resampled)'];
        
else
    csnew   = target;
    Rnew    = maprefcells(DEM.georef.XWorldLimits,DEM.georef.YWorldLimits,csnew,csnew);
    xlimitsn = Rnew.XWorldLimits;
    ylimitsn = Rnew.YWorldLimits;

    Z = imtransform(DEM.Z,T,method,...
        'Udata',[u(1) u(end)],'Vdata',[v(1) v(end)],...
        'Xdata',xlimitsn,'Ydata',ylimitsn,...
        'XYscale',[csnew csnew],...
        'FillValues',fillval);

    DEMr = GRIDobj(Z,Rnew);
    DEMr.georef.ProjectedCRS = DEM.georef.ProjectedCRS;

end

DEMr.name    = [DEM.name ' (resampled)'];

end
%% ----------------------------------------------------------------------
function mustBeGRIDobjOrPositiveScalar(x)

tf = isa(x,'GRIDobj') || (isscalar(x) && x > 0);
if ~tf 
    error("TopoToolbox:validation","Input must be GRIDobj or a positive scalar.")
end
end

%% ----------------------------------------------------------------------
% functions to create imref2d objects
function [R,Z] = GRIDobj2imref2d(DEM,isproj)
if isproj
    R = imref2d(DEM.georef.RasterSize,DEM.georef.XWorldLimits,DEM.georef.YWorldLimits);
    columnsStartNorth = strcmp(DEM.georef.ColumnsStartFrom,'north');
else
    if ~isempty(DEM.georef)
        R = imref2d(DEM.georef.RasterSize,DEM.georef.LongitudeLimits,DEM.georef.LatitudeLimits);
        columnsStartNorth = strcmp(DEM.georef.ColumnsStartFrom,'north');
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



