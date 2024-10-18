function I = isnan(DEM)

%ISNAN Returns array elements that are NaNs as logical grid
%
% Syntax
%
%     I = isnan(DEM)
%
% Description
%
%     overloaded isnan for GRIDobj. 
%
% Input arguments
%
%     DEM   GRIDobj
%
% Output arguments
%
%     I     GRIDobj with underlying type logical where true elements
%           indicate NaNs in the DEM
%
% See also: isnan
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 18. October, 2024

I = DEM;
I.Z = isnan(DEM.Z);
I.name = [DEM.name ' (isnan)'];