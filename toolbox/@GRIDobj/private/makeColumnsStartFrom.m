function DEM = makeColumnsStartFrom(DEM,direction)
%Forces column direction in GRIDobj to follow specified direction
%
% Syntax
%
%     DEM = makeColumnsStartFrom(DEM,'north')
%
% See also: GRIDobj, GRIDobj/reproject2utm, egm96heights, imtransform
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 31. March, 2025

arguments
    DEM GRIDobj
    direction {mustBeMember(direction,{'south','north'})} = 'north'
end

actualDirection = DEM.georef.ColumnsStartFrom;
if isequal(actualDirection,direction)
    return
end

R = DEM.georef;
R.ColumnsStartFrom = direction;
Z = flipud(DEM.Z);
DEM = GRIDobj(Z,R);



