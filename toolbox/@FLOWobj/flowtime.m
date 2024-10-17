function [OUT1, OUT2] = flowtime(FD,V,IX) 

%FLOWTIME Flow time (and distance) in upstream direction
%
% Syntax
%
%     [TIME, DIST] = flowtime(FD,V,IX)
%
% Description
%
%     flowtime calculates the horizontal distance D along the flow network
%     in upstream direction from a location IX, similar to flowdistance.m .
%     Moreover, it computes the flow time required to reach the location IX
%     based on a grid of flow velocity V.
%
% Input arguments
%
%     FD          flow direction (FLOWobj)
%     V           grid of flow velocity in m/s (GRIDobj)
%     IX          linear index of the seed
%
% Output arguments
%
%     T            time grid in s (GRIDobj)
%     D            distance grid (GRIDobj)
%
% See also: FLOWobj, FLOWobj/flowacc, GRIDobj
% 
% Author: Philippe Steer (philippe.steer[at]univ-rennes.fr)
% Based on flowdistance.m written by Wolfgang Schwanghart
% Date: 07. October, 2024

%% check input arguments

if ~strcmpi(FD.type,'single')
    error('TopoToolbox:FLOWobj','flowtime requires flow distance type to be single')
end

SEED = false(FD.size);
SEED(IX) = true;

%% Do calculation
ixtemp  = FD.ix;
ixctemp = FD.ixc;

DIST = getdistance(ixtemp,ixctemp,FD.size,FD.cellsize,'single');
TIME = DIST./V.Z(ixtemp);

cl   = class(DIST);

%% Upstream distance calculation
D     = inf(FD.size,cl);
D(SEED) = 0;
T     = inf(FD.size,cl);
T(SEED) = 0;
start = find(SEED(ixctemp),1,'last');
for r = start:-1:1
    D(ixtemp(r)) = min(D(ixctemp(r))+DIST(r),D(ixtemp(r)));
    T(ixtemp(r)) = min(T(ixctemp(r))+TIME(r),T(ixtemp(r)));
end
D(isinf(D)) = nan;
T(isinf(D)) = nan;        
%% Prepare Output
OUT1 = GRIDobj(FD); OUT2 = GRIDobj(FD);
% write time output to GRIDobj
OUT1.Z = T;
OUT1.zunit = '';
OUT1.name  = ['upstream flow time'];
% write distance output to GRIDobj
OUT2.Z = D;
OUT2.zunit = '';
OUT2.name  = ['upstream flow distance'];

