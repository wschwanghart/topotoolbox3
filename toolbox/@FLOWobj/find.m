function [ix,ixc,frac] = find(FD)

%FIND Find indices and values of edges in the flow direction graph
%
% Syntax
%
%     [ix,ixc] = find(FD); % for single flow direction
%     [ix,ixc,frac] = find(FD); % for multiple flow direction
%
% Description
%
%     FIND returns the linear indices of nodes in the flow network stored
%     in an instance FD of FLOWobj. frac contains the fraction that each
%     cell in ix delivers to its downstream neighbor in ixc. frac is an
%     empty array if FD represents single flow directions.
%
% Input arguments
%
%     FD     FLOWobj
%
% Output arguments
%
%     ix     giver indices
%     ixc    receiver indices
%     frac   fraction
%
% See also: FLOWobj, FLOWobj2M
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 24. June, 2024

arguments
    FD   FLOWobj
end

ix = FD.ix;
ixc = FD.ixc;
frac = FD.fraction;