function [DEMr,zone] = reproject2utm(DEM,res,varargin)

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
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 31. May, 2024



if ~isa(res,'GRIDobj')
    % get latitude and longitude vectors
    [lon,lat] = getcoordinates(DEM);
    % and calculate centroid of DEM. The centroid is used to
    % get the utmzone
    lonc = sum(lon([1 end]))/2;
    latc = sum(lat([1 end]))/2;
    zone        = utmzone(latc,lonc);
    if double(upper(zone(end)))>=78
        hemisphere = 'N';
    else
        hemisphere = 'S';
    end
    zone = str2double(zone(1:end-1));
else
    zone = [];
    hemisphere = [];
    
end

% parse input arguments 
p = inputParser;
validmethods = {'bicubic','bilinear','nearest','linear'}; 
p.FunctionName = 'GRIDobj/reproject2UTM';
% required
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
addRequired(p,'res',@(x) (~isa(x,'GRIDobj') && isscalar(x) && x > 0) || isa(x,'GRIDobj'));
% optional
addParameter(p,'zone',zone,@(x) ischar(x) || validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',60}));
addParameter(p,'hemisphere',hemisphere,@(x) ischar(validatestring(x,{'N','S'})))
addParameter(p,'method','linear',@(x) ischar(validatestring(x,validmethods)));

parse(p,DEM,res,varargin{:});

if isa(p.Results.res,'GRIDobj')
    DEMr = project(DEM,p.Results.res,'method',p.Results.method);
    zone = [];
else

    % Calculate EPSG-code from zone
    % If on southern hemisphere, 32700 + zone
    % If on northern hemisphere, 32600 + zone
    zone = p.Results.zone;
    if isstring(zone) | ischar(zone)
        zone = str2double(zone(1:end-1));

        if double(upper(zone(end)))>=78
            epsg = 32600 + zone;
        else
            epsg = 32700 + zone;
        end
    else
        switch lower(p.Results.hemisphere)
            case 'n'
                epsg = 32600 + zone;
            case 's'
                epsg = 32700 + zone;
        end
    end
    DEMr = project(DEM,epsg,'res',p.Results.res,'method',p.Results.method);
end

