function DEM = xyz2GRIDobj(xyz)

%XYZ2GRIDobj Convert xyz-array to GRIDobj
%
% Syntax
%
%     DEM = xyz2GRIDobj(xyz)
%
% Description
%
%     xyz2GRIDobj converts a set of scattered or gridded XYZ point data
%     into a raster-based GRIDobj, assuming that the x–y coordinates lie on
%     a regular grid.
%
% Input arguments
%
%     xyz     n × 3 numeric matrix with x-, y-coordinates and z-values
%
% Output arguments
%
%     DEM     GRIDobj
%
% Example
%
%     [X,Y] = meshgrid(1:25);
%     X = X(:);
%     Y = Y(:);
%     Z = X + Y + rand(size(X));
%     DEM = xyz2GRIDobj([X Y Z]);
%     imagesc(DEM)
%    
% See also: GRIDobj
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 6. June, 2025


arguments
    xyz (:,3) {mustBeNumeric}
end

x   = xyz(:,1);
y   = xyz(:,2);
z   = xyz(:,3);
z   = single(z);

[xi,~,locx]   = unique(x,'sorted');
[yi,~,locy]   = unique(y,'sorted');

dx = diff(xi);
dy = diff(yi);

tol = 1e-10;

isRegularX = max(abs(dx - dx(1))) < tol;
isRegularY = max(abs(dy - dy(1))) < tol;

if ~(isRegularY) || ~(isRegularX)
    error('The coordinates must be on a regular grid.')
end

Z  = nan(numel(yi),numel(xi),'single');
ix = sub2ind([numel(yi),numel(xi)],locy,locx);
Z(ix) = z;

DEM = GRIDobj(xi,yi,Z);