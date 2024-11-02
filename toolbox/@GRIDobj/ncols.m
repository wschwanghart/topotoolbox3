function nc = ncols(DEM)
%NCOLS Return the number of columns in a DEM
%
% Syntax
%   
%     nr = ncols(DEM)
%
% Description
%
%     The function returns the number of rows in a GRIDobj DEM.
%
% See also: GRIDobj/nrows
% 
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 2. November, 2024

nc = DEM.size(2);