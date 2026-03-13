function [P,a,locb] = removeduplicates(P)

%REMOVEDUPLICATES Remove duplicate points
%
% Syntax
%
%     P2 = removeduplicates(P)
%     [P2,a,locb] = ...
%
% Description
%
%     removeduplicates removes duplicate points in the point pattern P by
%     retaining only one of the duplicates.
%
% Input arguments
%
%     P      point pattern on stream network (class PPS)
%
% Output arguments
%
%     P2     point pattern on stream network (class PPS)
%     a      linear index into points, so that P2.PP = P.PP(a)
%     locb   linear index, so that P.PP = P2.PP(locb).
%
% See also: PPS, PPS/hasduplicates, unique
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 10. March, 2026

arguments
    P   PPS
end

[u,a,b] = unique(P.PP,'stable');

if numel(u) == numel(P.PP)
    return
end

if nargout >= 2
    ix = (1:npoints(P))';
    ixp = accumarray(b,ix,[numel(u) 1],@(x) {x});
    ixp = cellfun(@(x) sort(x),ixp,'UniformOutput',false);
    locb = cellfun(@(x) x(2:end),ixp,'UniformOutput',false);
end

P.PP = P.PP(a);


