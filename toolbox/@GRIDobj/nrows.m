function nr = nrows(DEM)
%NROWS Return the number of rows in a DEM
%
% Syntax
%   
%     nr = nrows(DEM)
%
% Description
%
%     The function returns the number of rows in a GRIDobj DEM.
%
% See also: GRIDobj/ncols
% 
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 2. November, 2024

nr = DEM.size(1);