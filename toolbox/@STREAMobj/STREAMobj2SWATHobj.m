function [SW] = STREAMobj2SWATHobj(S,DEM,options)

%STREAMOBJ2SWATHOBJ Create swath profile (SWATHobj) from stream network
%
% Syntax
%
%     SW = STREAMobj2SWATHobj(S,DEM,'pn','pv',...)
%
% Description
%
%     STREAMobj2SWATHobj creates a swath profile along individual reaches
%     of a stream network. If the SWATHobj was created from a STREAMobj
%     with multiple channels, the resulting SWATHobj's will be stored in
%     cells.
% 
% Input arguments
%
%     S     STREAMobj
%     DEM   digital elevation model (DEM)
%     pn,pv parameter name value pairs as used in SWATHobj
%
% Output arguments
%
%     SW    SWATHobj or cell array of SWATHobjs
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(trunk(S));
%     SW = STREAMobj2SWATHobj(S,DEM);
%     plotdz(SW);
%
%
% See also: SWATHobj
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
%         Wolfgang Schwanghart (schwangh[@]uni-potsdam.de)
% Date: 13. September, 2024

arguments
    S    STREAMobj
    DEM  GRIDobj
    options.width (1,1) = 1e4
    options.dx (1,1) = DEM.cellsize
    options.dy (1,1) = DEM.cellsize
    options.gap (1,1) = 0
    options.keepdist (1,1) = true
    options.keeptrace (1,1) = false
    options.keepnodes (1,1) = false
    options.smooth (1,1) = false
    options.smoothmethod = 'sgolay'
    options.hillshade (1,1) = false
    options.smoothingfactor (1,1) = .25
end

% validate further
validatealignment(S,DEM);

% Obtain individual tributaries 
CS = STREAMobj2cell(S,'tributaries');
netd = S.distance; % Upstream distance for the entire network

% Preallocate output
CSW = cell(size(CS));

options = namedargs2cell(options);

for r = 1:numel(CS)
    [x,y,d] = STREAMobj2XY(CS{r},nal2nal(CS{r},S,netd));
    x = flipud(x(1:end-1));
    y = flipud(y(1:end-1));
    d = flipud(d(1:end-1));
    
    CSW{r} = SWATHobj(DEM,x,y,d,options{:});
end

if length(CSW)==1
    SW = CSW{1};
else
    SW = CSW;
end


