function [OUT] = swath2latlon(SW)
%SWATH2LATLON Convert spatial fields in SWATHobj to geographic coordinates
%
% Syntax
%
%     OUT = swath2latlon(SW)
%
% Description
%
%     SWATH2LATLON uses the mapping toolbox to convert all spatial fields
%     in the SWATHobj SW to geographic coordinates (latitude, longitude).
%     Requires Mapping toolbox.
%
% Input arguments
%
%     SW     instance of SWATHobj
%
% Output arguments
%
%     OUT    instance of SWATHobj
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     SW = SWATHobj(DEM)
%     SW2 = swath2latlon(SW)
%
%
% Author: Dirk Scherler (scherler[at]caltech.edu)
% Date: June, 2024

if ~isProjected(SW)
    error('TopoToolbox:wrongInput',...
        'SWATHobj has a geographic or no projected coordinate system');
end

CRS = parseCRS(SW);
    
[SW.xy0(:,2),SW.xy0(:,1)] = projinv(CRS,SW.xy0(:,1),SW.xy0(:,2));
[SW.xy(:,2),SW.xy(:,1)] = projinv(CRS,SW.xy(:,1),SW.xy(:,2));
[y,x] = projinv(CRS,SW.X(:),SW.Y(:));
SW.X  = reshape(x,size(SW.X));
SW.Y  = reshape(y,size(SW.Y));

OUT   = SW;
OUT.georef = [];
OUT.georef.GeographicCRS = geocrs(4326);
    
