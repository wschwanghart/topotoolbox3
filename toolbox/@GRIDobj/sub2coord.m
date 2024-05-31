function [x,y] = sub2coord(DEM,r,c)

%SUB2COORD convert subscripts to x and y coordinates
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

r = r(:);
c = c(:);

xy    = (DEM.wf*double([c-1 r-1 ones(numel(r),1)]'))';
x = xy(:,1);
y = xy(:,2);