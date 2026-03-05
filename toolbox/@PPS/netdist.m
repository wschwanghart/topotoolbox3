function d = netdist(P,options)

%NETDIST Shortest network distance
%
% Syntax
%
%     d = netdist(P)
%     d = netdist(P,pn,pv,...)
%
% Description
%
%     netdist computes the distance along a stream network. It calculates
%     the distance of each node in the network to the nearest point in the
%     PPS object P.
%
% Input arguments
%
%     P       PPS
% 
%     Parameter name/value pairs
%
%     'dir'     'both' (default) calculates the distance in upstream and
%               downstream direction. 'up' calculates distances only in
%               upstream direction and 'down' in downstream direction.
% 
%     'split'   false (default) or true. If true, netdist will do the
%               computation in parallel on each drainage basin. If the
%               Parallel Computing Toolbox is available, that may be 
%               faster. However, splitting the network in its connected
%               components produces some considerable computational
%               overhead.
%
%     'distance' By default, distance is S.distance, the distance from the 
%               outlet. Yet, any other continuous increasing function can
%               be used, e.g. chi (see chitransform).
%
% Output arguments
%
%     d     node-attribute list with distances. Nodes that cannot be
%           reached will have a distance of inf. 
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S   = STREAMobj(FD,'minarea',500);
%     S   = klargestconncomps(S);
%     P   = PPS(S,'runif',10,'z',DEM);
%     d   = netdist(P);
%     tiledlayout(2,2)
%     nexttile; plotc(P.S,d); box on
%     nexttile; plotdz(P,'z',d); ylabel('Distance from point');
%     d   = netdist(P,'dir','up');
%     nexttile; plotc(P.S,d); box on
%     nexttile; plotdz(P,'z',d); ylabel('Distance from point');
%
% See also: STREAMobj/netdist, PPS/voronoi
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 5. March, 2026

arguments
    P  PPS
    options.dir = 'both'
    options.split (1,1) = false
    options.distance = P.S.distance
end

if npoints(P) == 0
    d = inf(size(P.S.x));
    return
end

options = namedargs2cell(options);
d = netdist(P.S,P.S.IXgrid(P.PP),options{:});