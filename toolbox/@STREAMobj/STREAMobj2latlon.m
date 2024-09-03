function [lat,lon,varargout] = STREAMobj2latlon(S,varargin)

%STREAMOBJ2LATLON convert instance of STREAMobj to NaN-separated geographic coordinates
%
% Syntax
%
%     [lat,lon] = STREAMobj2latlon(S)
%     [lat,lon,a,b,...] = STREAMobj2latlon(S,A,B,...)
%
% Description
%
%     STREAMobj2latlon returns the NaN-punctuated latitude and longitude
%     vectors which can be used for easy plotting using the plot function.
%     With additional input arguments that are instances of GRIDobj,
%     additional vectors are generated with the respective values of the
%     grids A, B, etc. at the node locations of the stream network S.
%
%     Note that this function *requires the mapping toolbox*.
%
% Input arguments
%
%     S       streams (class STREAMobj)
%     A,B,... grids (class GRIDobj) or node attributes (e.g. as returned by
%             the function STREAMobj/streamorder
%
% Output arguments
%
%     lat,lon  coordinate vectors 
%     a,b,...  grid values at locations specified by x and y.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     [lat,lon] = STREAMobj2latlon(S);
%     plot(lon,lat)
%
% See also: STREAMobj/STREAMobj2XY
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 12. June, 2024

[CRS,isproj] = parseCRS(S);
if isnan(isproj)
    error("S has an undefined CRS.")
end

if ~isproj 
    error("S already has a geographic CRS.")
else
    
end
    
nr = numel(varargin);
c  = cell(1,nr);
[x,y,c{:}] = STREAMobj2XY(S,varargin{:});

[lat,lon] = projinv(CRS,x,y);
varargout = c;