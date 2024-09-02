function [x,y] = ind2coord(FD,ix)

%IND2COORD Convert linear index to x and y coordinates
%
% Syntax
%
%     [x,y] = ind2coord(DEM,ix)
%
% Description
%
%     ind2coord converts a linear index into an instance of GRIDobj to x
%     and y coordinates. The function returns column vectors.
%
% Input arguments
%
%     DEM     instance of GRIDobj
%     ix      linear index
%
% Output arguments
%
%     x,y     x- and y-coordinates  
%
% See also: GRIDobj/coord2ind, GRIDobj/getcoordinates
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 30. May, 2024


% Get intrinsic coordinates
[r,c] = ind2sub(FD.size,ix(:));
% Calculate map coordinates
xy    = (FD.wf*double([c-1 r-1 ones(numel(ix),1)]'))';
% Split coordinates in x and y vector
x = xy(:,1);
y = xy(:,2);