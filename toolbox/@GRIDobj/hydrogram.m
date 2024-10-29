function [OUT1,OUT2] = hydrogram(TIME,P,DT)

%HYDROGRAM Generate an hydrogram
%
% Syntax
%
%     [Q,t] = hydrogram(TIME,P,DT)
%
% Description
%
%     hydrogram computes the time evolution of discharge as an hydrogram
%     based on the flowtime TIME grid (see FLOWobj/flowtime.m) and the grid
%     of rainfall P. This function assumes that rainfall is a punctual
%     event corresponding to a dirac delta or impulse in time
%
% Input arguments
%
%     TIME  grid of flow time (GRIDobj)
%     P     grid of rainfall in m (GRIDobj)
%     DT    time interval of sampling (scalar) 
%
% Output arguments
%
%     Q      discharge array (temporal evolution)
%     t      time array
%
% Author: Philippe Steer (philippe.steer[at]univ-rennes.fr)
% Date: 07. October, 2024

%% check input arguments
% to be done

%% Compute results
% Find indices where T are defined
ind=find(TIME.Z>=0 & ~isinf(TIME.Z) & ~isnan(TIME.Z));
% Generate arrays of flow time and amplitude
timevec=TIME.Z(ind);
amplitudevec=P.Z(ind).*TIME.cellsize;
% Generate array of edges and bins for a temporal histogram
tedges=[0:DT:max(timevec)+DT];
tbins=(tedges(2:end)+tedges(1:end-1))./2;
% Compute an histogram of flow times
[~, ~, ibin] = histcounts(timevec, tedges);
% Accumulate discharge for each interval
Q = accumarray(ibin, amplitudevec);

%% Prepare Output
% Array of discharge
OUT1=Q;
% Array of time
OUT2=tbins;