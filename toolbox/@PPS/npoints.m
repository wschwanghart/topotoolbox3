function np = npoints(P)

%NPOINTS Number of points in the point pattern
%
% Syntax
%
%     np = npoints(P)
%
% Description
%
%     npoints returns the number of points in the the point pattern.
%
% Input arguments
%
%     P      point pattern on stream network (class PPS)
%
% Output arguments
%
%     np     number of points (integer scalar)
%
%
% See also: PPS, PPS/tlength 
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 12. November, 2024

arguments
    P PPS
end

np = numel(P.PP);