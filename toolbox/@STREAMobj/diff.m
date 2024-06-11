function dz = diff(S,z,fillval)

%DIFF Differences between adjacent pixels in a stream network
%
% Syntax 
% 
%     dz = diff(S,z)
%     dz = diff(S,z,fillval)
%
% Description
%
%     DIFF calculates the difference of each pixel i and its downstream
%     neighbor j so that dz(i) = dz(i)-dz(j). dz is a node-attribute list.
%     By default, stream nodes without downstream neighbor (outlets) receive a value of
%     0. If this value should be different, also set the input argument
%     fillval.
%
% Input arguments
%
%     S        STREAMobj
%     z        node-attribute list
%     fillval  value that dz should have at outlets (default = 0)
%
% Output arguments
%
%     dz    node-attribute list with differences
%
% Example: Calculate elevation offsets along stream network
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,'minarea',1000);
%     S = trunk(klargestconncomps(S));
%     dz = diff(S,DEM);
%     % You can calculate elevation using cumsum
%     zz = cumsum(S,dz,'upstream');
%     
% 
% See also: STREAMobj/gradient, STREAMobj/cumsum
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 11. June, 2024

arguments
    S   STREAMobj
    z   {mustBeGRIDobjOrNal(z,S)}
    fillval = 0
end

% get node attribute list with elevation values
z = ezgetnal(S,z);

dz = getnal(S) + fillval;
dz(S.ix) = z(S.ix)-z(S.ixc);