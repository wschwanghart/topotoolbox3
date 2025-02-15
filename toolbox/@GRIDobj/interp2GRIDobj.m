function A = interp2GRIDobj(DEM,x,y,v,method,extrapmethod)

%INTERP2GRIDobj Interpolate scattered data to GRIDobj
%
% Syntax
%
%     A = interp2GRIDobj(DEM,x,y,v)
%     A = interp2GRIDobj(DEM,x,y,v,method)
%     A = interp2GRIDobj(DEM,x,y,v,method,extrapolationmethod)
%
% Description
%
%     This function interpolates scattered data at the locations x,y with  
%     values v to a new GRIDobj spatially aligned with DEM. The function
%     wraps the built-in scatterInterpolant class and allows for the same
%     interpolation methods and extrapolation methods as this function.
%     Interpolation methods are 'linear', 'nearest', and 'natural'.
%     Extrapolation methods are 'linear', 'nearest' or 'none'.
%
% Input arguments
%
%     DEM      GRIDobj
%     x,y      coordinate vectors
%     v        vector with values
%     method   'linear' (default), 'nearest', or 'natural'
%     extrapolationmethod   'none' (default), 'linear', or 'nearest'
%
% Output argument
%
%     A        GRIDobj aligned with DEM with gridded values
%
% Example
%   
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     IX  = randperm(prod(DEM.size),200);
%     [x,y] = ind2coord(DEM,IX);
%     v = ((x-min(x))./1e5).^2 - ((y-min(y))./1e5).^1.5;
%     A = interp2GRIDobj(DEM,x,y,v,'natural','linear');
%     imageschs(DEM,A)
%     hold on
%     scatter(x,y,10,v,'filled','markeredgecolor','k')
%
% See also: scatteredInterpolant, GRIDobj/interp, GRIDobj
%        
% Author:  Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 8. October, 2024 

arguments
    DEM   GRIDobj
    x     {isnumeric}
    y     {isnumeric}
    v
    method {mustBeMember(method,{'linear','nearest','natural'})} = 'linear'
    extrapmethod ...
        {mustBeMember(extrapmethod,{'linear','nearest','none'})} = 'none'
end

if ~isequal(size(x),size(y),size(v))
    error('TopoToolbox:interp2GRIDobj','The input vectors must have same size')
end

% get coordinates from DEM
[xi,yi] = getcoordinates(DEM);

% create copy of DEM
A   = DEM;

% interpolate using scatteredInterpolant
F   = scatteredInterpolant(x,y,v,method,extrapmethod);
Z   = F({xi,yi});

% Transpose output
Z   = Z';

% Create output
A.Z = Z;
A.name = 'interp2GRIDobj';

