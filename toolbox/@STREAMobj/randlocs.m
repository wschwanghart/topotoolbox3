function varargout = randlocs(S,n,wr)

%RANDLOCS Random locations along the stream network
%
% Syntax
%
%     IX = randlocs(S,n)
%     IX = randlocs(S,n,wr)
%     [x,y] = randlocs(...)
%
% Descriptions
%
%     randlocs returns n random locations on the stream network S. You can
%     choose between sampling with replacement (wr = true) or without
%     replacement (wr = false, default).
%
% Input arguments
%
%     S      STREAMobj
%     n      number of locations. Note that when sampling without replacement
%            n may not be greater than the number of nodes in the network.
%     wr     true: with replacement
%            false (default): without replacement
%
% Output arguments
%
%     IX     linear index into the DEM
%     x,y    coordinates
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S  = STREAMobj(FD,'minarea',1e6,'unit','map');
%     [x,y] = randlocs(S,100);
%     plot(S)
%     hold on
%     plot(x,y,'s','MarkerFaceColor',[.6 .6 .6])
%
% See also: STREAMobj, PPS
% 
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 3. July, 2024

arguments
    S    STREAMobj
    n    (1,1) {mustBeInteger,mustBePositive} = 100
    wr   (1,1) = false
end

if wr
    idx = randi(numel(S.x),n,1);
else
    idx = randperm(numel(S.x),n);    
end

if nargout == 1
    varargout{1} = S.IXgrid(idx);
elseif nargout == 2
    varargout{1} = S.x(idx);
    varargout{2} = S.y(idx);
end


