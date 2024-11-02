function [x,y] = sub2coord(DEM,r,c)

%SUB2COORD Convert subscripts to x and y coordinates
%
% Syntax
%
%     [x,y] = sub2coord(DEM,r,c)
%
% Description
%
%     sub2coord converts a subscripts into an instance of GRIDobj to x
%     and y coordinates.
%
% Input arguments
%
%     DEM     instance of GRIDobj
%     r,c     row and column indices
%
% Output arguments
%
%     x,y     x- and y-coordinates (column vectors)
%
% See also: GRIDobj/coord2sub, GRIDobj/getcoordinates
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 30. May, 2024

arguments
    DEM  GRIDobj
    r
    c 
end

if numel(r) ~= numel(c)
    error('TopoToolbox:wronginput','r and c must have the same number of elements')
end

mustBeInRange(r,1,nrows(DEM))
mustBeInRange(c,1,ncols(DEM))

r = r(:);
c = c(:);

xy    = (DEM.wf*double([c-1 r-1 ones(numel(r),1)]'))';
x = xy(:,1);
y = xy(:,2);
