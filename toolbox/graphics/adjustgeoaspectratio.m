function adjustgeoaspectratio(h)
%ADJUSTGEOASPECTRATIO Adjusts data aspect ratio for latitude-longitude axes
%
% Syntax
%
%     adjustgeoaspectratio(h)
%
% Description
%
%     This function modifies the data aspect ratio of the specified axes to
%     reduce distortions in metric distances when the axes represent
%     latitude and longitude coordinates. It accounts for variations in
%     distance due to latitude, ensuring a more accurate spatial
%     representation.
%    
% Input arguments
%
%     h    axes handle
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 27. February, 2025

arguments
    h = gca
end

% Get axes limits
xlim = h.XLim;
ylim = h.YLim;
% Centers of both axes
mx    = sum(xlim)/2;
my    = sum(ylim)/2;
% Metric distances
dx   = distance(my,xlim(1),my,xlim(2),wgs84Ellipsoid);
dy   = distance(ylim(1),mx,ylim(2),mx,wgs84Ellipsoid);
% Metric distance per longitude/latitude degree
dx   = dx/range(xlim);
dy   = dy/range(ylim);
% Ratio
r    = dx/dy;
% Adjust data aspect ratio
h.DataAspectRatio = [1 r 1];