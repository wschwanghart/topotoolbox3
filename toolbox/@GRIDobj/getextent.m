function ext = getextent(DEM,latlonout)

%GETEXTENT return extent of a GRIDobj
%
% Syntax
%
%     ext = getextent(DEM,latlonout)
%
% Description
%
%     getextent returns the horizontal extent of the GRIDobj DEM. The
%     function will return the maximum extent in WGS84 geographical
%     coordinates if latlonout is true.
%
% Input arguments
%
%     DEM        GRIDobj
%     latlonout  false (default) or true. If true, getextent will return
%                the extent in geographical coordinates. This is, however,
%                only possible if DEM has a known projected coordinate 
%                system (DEM.georef.ProjectedCRS must be set), and if the
%                mapping toolbox is available. 
%
% Output arguments
%
%     ext        four element row vector with following format 
%                [west east south north]
%
%
% See also: GRIDobj, readopentopo, GRIDobj/getoutline
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 9. October, 2018

arguments
    DEM GRIDobj
    latlonout (1,1) {mustBeInteger} = 0
end

% Get the extent of the projected coordinates
[x,y] = getcoordinates(DEM);
csh   = DEM.cellsize/2;
ext   = [min(x)-csh max(x)+csh min(y)-csh max(y)+csh];

if latlonout
    % This requires the mapping toolbox and that the DEM has a projcrs
    [lat,lon] = projinv(DEM.georef.ProjectedCRS,...
                            [ext(1) ext(1) ext(2) ext(2)]',...
                            [ext(3) ext(4) ext(3) ext(4)]');
    ext   = [min(lon) max(lon) min(lat) max(lat)];
end
    
    

