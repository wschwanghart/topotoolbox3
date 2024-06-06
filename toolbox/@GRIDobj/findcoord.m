function varargout = findcoord(DEM)

%FINDCOORD Find coordinates of nonzero elements in GRIDobj
%
% Syntax
%
%     [x,y] = findcoord(DEM)
%     [lon,lat] = findcooord(DEM)
%     [...,val] = findcoord(DEM)
%
% Description
%
%     This function returns the coordinates of nonzero values in DEM. If 
%     the DEM is in a geographic coordinate system, then the function
%     returns lon and lat.
%
% Input arguments
%
%     DEM     GRIDobj
%
% Output arguments
%
%     x,y     coordinate pair
%     val     value of nonzero element in DEM.Z
%
% See also: find, GRIDobj/coord2ind, GRIDobj/ind2coord, GRIDobj/find           
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 4. June, 2024

ix = find(DEM);
[x,y] = ind2coord(DEM,ix);

varargout{1} = x;
varargout{2} = y;
if nargout == 3
    varargout{3} = DEM.Z(ix);
end

