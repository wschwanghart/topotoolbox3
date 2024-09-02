function tf = ismulti(FD,typetest)
%ISMULTI Determine whether FD is multi or single flow direction
%
% Syntax
%
%     tf = ismulti(FD)
%
% Description
%
%     ISMULTI returns true if a FLOWobj contains multiple flow directions
%
% Input arguments
%
%     FD     FLOWobj
%
% Output arguments
%
%     tf     true or false
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     ismulti(FD)
% 
%     ans =
% 
%       logical
% 
%        0
%
% See also: FLOWobj
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 31. August, 2024

if nargin == 1
    typetest = true;
end

if typetest
    tf = strcmp(FD.type,'multi');
else 
    tf = any(histcounts(FD.ix,1:(prod(FD.size)+1))>1);
end
