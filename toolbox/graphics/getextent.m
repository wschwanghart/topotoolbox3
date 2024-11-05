function extent = getextent(ax)

%GETEXTENT Get current axis extent 
%
% Syntax
%
%     extent = getextent(ax)
%
% Description
%
%     getextent and setextent are two small wrapper functions (around set
%     and get) to quickly zoom to a specified extent in an axis object.
%
% Example
% 
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     imageschs(DEM);
%     % execute til here and zoom to a desired level using the zoom tool
%     e = getextent;
%     imageschs(DEM,gradient8(DEM))
%     setextent(e)
%
% See also: IMAGESCHS, SETEXTENT, PADEXTENT
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 4. November 2024

arguments
    ax {mustBeA(ax,'matlab.graphics.axis.Axes')} = gca
end

extent = get(ax,{'xlim','ylim'});  % Get axes limits.

