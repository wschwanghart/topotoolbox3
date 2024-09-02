function stats = drainagebasinstats(FD,L,varargin)

%DRAINAGEBASINSTATS Zonal statistics on drainage basins
% 
% Syntax
%
%     stats = drainagebasinstats(FD,IX)
%     stats = drainagebasinstats(FD,L)
%     stats = drainagebasinstats(FD,MS)
%     stats = drainagebasinstats(...,'name1',GRIDobj1,...)
%
% Description
%
%     drainagebasinstats calculates zonal statistics of one or several 
%     grids based on drainage basins derived from a FLOWobj. It differs
%     from using drainage basin grids or other label matrices by allowing
%     drainage basins to be nested. 
%
% Input arguments
%
%     FD         FLOWobj
%     IX         linear index into the DEM (GRIDobj) from which FD was
%                derived.
%     L          label grid (GRIDobj). 
%     MS         polygon mapping structure array
%
%     Additional parameter name/value pairs
%
%     'waitbar'               {false} or true
%     'geotable'              {true} or false. If true, the function
%                             returns a geotable. Otherwise, a mapping 
%                             structure will be returned.
%     'name1',GRIDobj1,...    Concatenate pairs of grid names and GRIDobj
%                             that shall be used to calculate statistics.
%                             'name1' will be used as prefix in the output
%                             structure array stats. Be sure to use unique
%                             names.
%
% Output arguments
%
%     stats     structure array with statistics
%
% Example: Plot the distribution of mean gradients upstream of channelheads
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     S = STREAMobj(FD,'minarea',1000);
%     ix = streampoi(S,'channelheads','ix');
%     G  = gradient8(DEM);
%     stats = drainagebasinstats(FD,ix,'grad',G);
%     histogram(stats.grad_mean)
%
% See also: FLOWobj, drainagebasins
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 31. August, 2024

p = inputParser;
p.FunctionName = 'FLOWobj/drainagebasinstats';
p.KeepUnmatched = true;

addParameter(p,'waitbar',false)
addParameter(p,'geotable',true)
parse(p,varargin{:});

% check additional input grids
C = namedargs2cell(p.Unmatched);
narg = numel(C)/2;


if isa(L,'GRIDobj')
    % get unique labels in L
    l = unique(L.Z(:));

    ii = ~isnan(l) & l~=0;
    l = l(ii);
    nrlabels = numel(l);
    inp = 'GRIDobj';
elseif isstruct(L)
    inp = 'shape';
    nrlabels = numel(L);
    l   = 1:nrlabels;
    stats = L;
else
    inp = 'index';
    ix  = L;
    nrlabels = numel(ix);
    l   = 1:nrlabels;
end

wb = nrlabels > 2 && p.Results.waitbar;
if wb
    h = waitbar(0,'please wait');
end

for r=1:nrlabels
    if wb
        waitbar(r/nrlabels,h);
    end
    switch inp
        case 'GRIDobj'
            I = L==l(r);
        case 'index'
            I = GRIDobj(FD,'logical');
            I.Z(ix(r)) = true;
        case 'shape'
            I = GRIDobj(FD,'logical');
            RR = I.wf';
            P = [L(r).X(:) - I.wf(1,3), L(r).Y(:) - I.wf(2,3)]/RR(1:2,:);
            x = P(:,1);
            y = P(:,2);
            
            ii = isnan(x) | isnan(y);
            x = x(~ii);
            y = y(~ii);
            
            minx = floor(min(x));
            maxx = ceil(max(x));
            miny = floor(min(y));
            maxy = ceil(max(y));
            
            ii = poly2mask(x-minx+1,y-miny+1,maxx-minx+1,maxy-miny+1);
            
            if any(ii)
                I.Z(minx:maxx,miny:maxy) = ii;
            else
                IX = coord2ind(I,L(r).X(:),L(r).Y(:));
                IX(isnan(IX)) = [];
                IX = unique(IX);
                I.Z(IX) = true;
            end
    end
    
    I = dependencemap(FD,I);% & ~I;
    
    % In addition, optionally get outline
    stats(r).Geometry = 'Polygon';
    [~,X,Y] = GRIDobj2polygon(I);
    stats(r).X = X;
    stats(r).Y = Y;
    
    % Some basic geometric values
    stats(r).label = double(l(r));
    stats(r).upslopearea  = nnz(I.Z)*I.cellsize^2;
    [xx,yy] = findcoord(I);
    stats(r).Xcentr = mean(xx);
    stats(r).Ycentr = mean(yy);
    
    if narg > 0
        for r2 = 1:narg
            v = C{(r2-1)*2 +2}.Z(I.Z);
            v = v(:);
            v = v(~isnan(v));
            
            if islogical(v)
                stats(r).([C{(r2-1)*2 +1} '_mean']) = mean(double(v));
            else
            v = double(v);
            stats(r).([C{(r2-1)*2 +1} '_mean']) = mean(v);
            stats(r).([C{(r2-1)*2 +1} '_std']) = std(v);
            prc = prctile(v,[1 5 33 50 66 95 99]);
            stats(r).([C{(r2-1)*2 +1} '_prc01']) = prc(1);
            stats(r).([C{(r2-1)*2 +1} '_prc05']) = prc(2);
            stats(r).([C{(r2-1)*2 +1} '_prc33']) = prc(3);
            stats(r).([C{(r2-1)*2 +1} '_median']) = prc(4);
            stats(r).([C{(r2-1)*2 +1} '_prc66']) = prc(5);
            stats(r).([C{(r2-1)*2 +1} '_prc95']) = prc(6);
            stats(r).([C{(r2-1)*2 +1} '_prc99']) = prc(7);
            
            minv = min(v); if isempty(minv); minv = nan; end
            maxv = max(v); if isempty(maxv); maxv = nan; end
            stats(r).([C{(r2-1)*2 +1} '_min']) = minv;
            stats(r).([C{(r2-1)*2 +1} '_max']) = maxv;
            stats(r).([C{(r2-1)*2 +1} '_range']) = maxv-minv;
            stats(r).([C{(r2-1)*2 +1} '_kurtosis']) = kurtosis(v);
            stats(r).([C{(r2-1)*2 +1} '_skewness']) = skewness(v);
           
            
            end
        end
    end
end

% Save as geotable, if required.
if p.Results.geotable
    stats = mapstruct2geotable(stats,'coordinateSystemType','planar',...
        'CoordinateReferenceSystem',FD.georef.ProjectedCRS);
end

if wb
    close(h);
end

