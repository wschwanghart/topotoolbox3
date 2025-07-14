function varargout = hypscurve(DEM,bins,options)

%HYPSCURVE plot hypsometric curve of a digital elevation model
%
% Syntax
%
%     hypscurve(DEM)
%     hypscurve(DEM,nrbins)
%     ax = hypscurve(...)
%     [rf,elev] = hypscurve(...)
%
% Description
%
%     A hypsometric curve is an empirical cumulative distribution function
%     of elevations in a catchment. hypscurv plots the hypsometric curve or
%     returns the relative frequencies of elevation values. Optionally, 
%     hypscurve bins the data in equally spaced containers. 
%
% Input
%
%     DEM       digital elevation model (GRIDobj)
%     nrbins    number of bins
%
%     Parameter name/value pairs
%
%     'clip'    STREAMobj to return only the hypsometrical curve of the
%               stream network
%
% Output
%
%     ax        axis handle
%     rf        relative frequencies (in %)
%     elev      elevation
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     hypscurve(DEM,50)
%
% See also: histogram
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 11. July, 2025

arguments
    DEM    GRIDobj
    bins = []
    options.clip = []
end

if isempty(options.clip)
    z = DEM.Z;
    z = z(~isnan(z));
    z = z(:);
else
    if isa(options.clip,'STREAMobj')
        z = getnal(options.clip,DEM);
        z = z(~isnan(z));
    end
end

if ~isempty(bins)
    [n,edges] = histcounts(z,bins);
    elev = edges(1:end-1) + diff(edges)/2;
    n = flipud(n(:));
    elev = flipud(elev(:));
    n = cumsum(n);
    linestyle = '-';
else
    elev = sort(z,'descend');
    
    % return unique entries only
    I = [true;diff(elev)~=0] | [flipud(diff(flipud(elev))~=0); true];
    elev = elev(I);
    n    = find(I);
    linestyle = '-';
 
end

% get relative frequencies
n = n./n(end) * 100;

% plot results
if nargout ~= 2
    axis_handle = plot(n,elev,linestyle);
    axis xy
    xlabel('Compl. cum. frequency [%]')
    ylabel('Elevation [m]')
end

% prepare output
if nargout == 1
    varargout{1} = axis_handle;
elseif nargout == 2
    varargout{1} = n;
    varargout{2} = elev;
end
    



