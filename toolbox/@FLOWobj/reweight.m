function FD = reweight(FD,DEM,beta,options)

%REWEIGHT Re-calculate flow weights in a multiple flow direction FLOWobj
%
% Syntax
%
%     FD = reweight(FD,DEM,beta)
%     FD = reweight(FD,DEM,beta,pn,pv,...)
%
% Description
%
%     The function reweight takes a multiple flow direction object FD and
%     re-calculates the proportions at which flow is distributed from each  
%     pixel to its downstream neighbors. By default, proportions in FD 
%     calculated with FLOWobj(DEM,'multi') are proportional to slope, so 
%     that beta = 1.
%
%     Reweight recalculates the slopes based on the DEM and then uses beta
%     to reweight and normalizes the edges of the flow network. This is the
%     approach taken by Holmgren (1994).
%
%     Note that reweight does not change the topology of the flow network
%     so that no links are added nor removed. If slopes between two nodes
%     in the network are zero or negative, their slope is adjusted to have
%     a very small value (0.1 * the smallest non-negative slope value).
%
% Input arguments
%
%     FD      FLOWobj (multiple flow directions)
%     DEM     Digital elevation model (GRIDobj)
%     beta    exponent >= 0 
%
%     Parameter name/value pairs
%
%     minslope    scalar (0.1). This indicates the multiplier used to 
%                 derive the minimum slope value for links that have zero
%                 or negative slopes.
%     stable      false or true. Uses a numerically stable approach when
%                 using large exponents. By default false, but if beta > 5,
%                 then the stable approach is chosen.
%     
% Output arguments
%
%     FD      FLOWobj
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'multi');
%     FD1 = reweight(FD,DEM,0);
%     FD2 = reweight(FD,DEM,10);
%     A1 = flowacc(FD1);
%     A2 = flowacc(FD2);  
%     tiledlayout(1,2,"TileSpacing","compact");
%     ax(1) = nexttile; imageschs(DEM,log(A1),'colormap',flowcolor);
%     ax(2) = nexttile; imageschs(DEM,log(A2),'colormap',flowcolor);
%     linkaxes(ax,'xy');
%     ext = {[3.7769e+05 3.8353e+05]   [3.7911e+06 3.7948e+06]};
%     setextent(ext,ax(1))
%
% See also: FLOWobj, FLOWobj/multi_normalize
%  
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 6. April, 2026

arguments 
    FD FLOWobj
    DEM GRIDobj
    beta (1,1) {mustBeNonnegative}
    options.minslope (1,1) {mustBeNonnegative} = 0.1
    options.stable (1,1) = beta > 5
end

[X,Y] = getcoordinates(DEM,'mat');
dz = double(DEM.Z(FD.ix)) - double(DEM.Z(FD.ixc));
dx = sqrt((X(FD.ix) - X(FD.ixc)).^2 + ...
          (Y(FD.ix) - Y(FD.ixc)).^2);
FD.fraction = dz./dx;
minNonZeroFrac = min(FD.fraction(FD.fraction>0));
FD.fraction = max(FD.fraction,minNonZeroFrac*0.1);

if ~options.stable
    FD.fraction = FD.fraction.^beta;
    FD = multi_normalize(FD);
else
    % [~,~,locb] = unique(FD.ix);
    % [W,maxY] = splitapply(@(x) normfun(x,beta),FD.fraction,locb);
    % FD.fraction = exp(beta.*log(FD.fraction) - maxY(locb))./W(locb);
    W = accumarray(FD.ix,FD.fraction,[numel(DEM.Z) 1],@(x) normfun(x,beta));
    maxY = accumarray(FD.ix,FD.fraction,[numel(DEM.Z) 1],@(x) normfun2(x,beta));    
    FD.fraction = exp(beta.*log(FD.fraction) - maxY(FD.ix))./W(FD.ix);

end

end

function [sumw,maxy] = normfun(x,beta)
    y = beta * log(x);
    maxy = max(y);
    y = y - maxy;
    w = exp(y);
    % w = w / sum(w);
    sumw = sum(w);
end

function maxy = normfun2(x,beta)
    y = beta * log(x);
    maxy = max(y);
end