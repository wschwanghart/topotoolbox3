function [x,y] = getcoordinates(DEM,type)

%GETCOORDINATES get coordinate vectors of an instance of GRIDobj
%
% Syntax
%
%     [x,y] = getcoordinates(DEM)
%     [x,y] = getcoordinates(DEM,type)
%
% Input arguments
%
%     DEM    grid (class: GRIDobj)
%     type   'vector' (default). Alternatively, you can return coordinate
%            matrices ('matrix') or GRIDobjs ('GRIDobj')
%
% Output arguments
%
%     x      coordinate vector in x direction (row vector)
%     y      coordinate vector in y direction (column vector)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     [x,y] = getcoordinates(DEM);
%     surf(x,y,double(DEM.Z))
%     axis image; shading interp; camlight
%     
%
%
% See also: GRIDobj2mat
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 27. November, 2025

[x,y] = wf2XY(DEM.wf,DEM.size);

if nargin == 1
    type = 'vector';
else
    type = validatestring(type,{'vector','matrix','GRIDobj'});
end

switch type
    case 'vector'
        x = x(:)';
        y = y(:);
    case 'matrix'
        [x,y] = meshgrid(x,y);
    case 'GRIDobj'
        [x,y] = meshgrid(x,y);
        X = GRIDobj(DEM);
        X.Z = x;
        Y = GRIDobj(DEM);
        Y.Z = y;
        x = X;
        y = Y;
end