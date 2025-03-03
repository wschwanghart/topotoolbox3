function [DEMr,zone] = reproject2utm(DEM,res,options)

%REPROJECT2UTM Reproject DEM with WGS84 coordinate system to UTM-WGS84 
%
% Syntax
% 
%     [GRIDr,zone] = reproject2utm(GRID,res)
%     [GRIDr,zone] = reproject2utm(GRID,res,pn,pv,...)
%     GRIDr        = reproject2utm(GRID,GRID2)
%     GRIDr        = reproject2utm(GRID,GRID2,'method',method)
%
% Description
%
%     Reproject a grid (GRIDobj) with WGS84 geographic coordinates to UTM 
%     WGS84 (requires the mapping toolbox and image processing toolbox).
%
% Input arguments
%
%     GRID     raster (GRIDobj) with WGS84 geographic coordinates. 
%     res      spatial resolution in x- and y-direction (scalar)
%     GRID2    raster (GRIDobj) with projected coordinate system to which
%              GRID shall be projected. The resulting grid will be
%              perfectly spatially aligned (same cellsize, same upper left
%              egde, same size) with GRID2.
%     
% Parameter name/value pairs
%
%     zone       is automatically determined. If supplied, the value must
%                be a string, e.g., '32T'. Note that this function requires
%                the full grid zone reference that includes the uppercase
%                letter indicating the latitudinal band. 
%     method     interpolation method ('bilinear' (default), 'bicubic',  
%                or 'nearest')
%
% Output arguments
%
%     GRIDr    raster (GRIDobj) with UTM-WGS84 projected coordinates
%     zone     utm zone (string)
%
%
% See also: GRIDobj, imtransform, utmzone, GRIDobj/project
%           
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 27. February, 2025

arguments
    DEM  GRIDobj
    res
    options.zone = []
    options.method {mustBeMember(options.method,{'linear','bilinear','bicubic','nearest'})} = 'bilinear'
    options.hemisphere {mustBeMember(options.hemisphere,{'N','S',''})} = ''
end

% If no resolution but a GRIDobj is provided, that's easy. Just invoke
% project
if isa(res,'GRIDobj')
    DEMr = project(DEM,res,'method',options.method);
    zone = [];
    return
end


% Otherwise, we need to determine the epsg-code of the zone 

% get latitude and longitude vectors
[lon,lat] = getcoordinates(DEM);

% and calculate centroid of DEM. The centroid is used to
% get the utmzone
lonc = sum(lon([1 end]))/2;
latc = sum(lat([1 end]))/2;
    
% If no zone is provided, then we'll determine it using utmzone
if isempty(options.zone)
    zone        = utmzone(latc,lonc);
 
    if double(upper(zone(end)))>=78
        hemisphere = 'N';
    else
        hemisphere = 'S';
    end
    zone = str2double(zone(1:end-1));
else
    zone = str2double(options.zone(1:end-1));
    if double(upper(options.zone(end)))>=78
        hemisphere = 'N';
    else
        hemisphere = 'S';
    end
end

% Any hemisphere preferred?
if ~isempty(options.hemisphere)
    hemisphere = options.hemisphere;
end

% Calculate epsg-code
switch lower(hemisphere)
    case 'n'
        epsg = 32600 + zone;
    case 's'
        epsg = 32700 + zone;
end

DEMr = project(DEM,epsg,'res',res,'method',options.method);


