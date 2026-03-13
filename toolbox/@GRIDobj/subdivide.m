function CDEM = subdivide(DEM,options)

%SUBDIVIDE Subdivide a GRIDobj into a cell array of GRIDobjs
%
% Syntax
%
%     CDEM = subdivide(DEM)
%     CDEM = subdivide(DEM,pn,pv,...)
%
% Description
% 
%     The function subdivides a GRIDobj DEM into tiles which are stored in
%     a cell array of GRIDobj. 
%     
% Input arguments
%
%     DEM      GRIDobj
%
%     Parameter name/value pairs
%     
%     tilesize  [100 100] (1x2) array of tile size in  
%
% Output argument
%
%     CDEM    Cell array of GRIDobj
%
% Example
%
%     DEM = GRIDobj("srtm_bigtujunga30m_utm11.tif");
%     CDEM = subdivide(DEM,"tilesize",[200 200]);
%     tiledlayout("flow"); 
%     for r = 1:numel(CDEM); 
%          nexttile; 
%          imageschs(CDEM{r},'colorbar',false,'ticklabels','none'); 
%     end
%     
%
% See also: GRIDobj/crop, GRIDobj/cropbyregion
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 2. March, 2026

arguments
    DEM
    options.tilesize = [100 100]
end

cs = DEM.cellsize;
L = reclabel(DEM,cs*options.tilesize(2),cs*options.tilesize(1));
CDEM = cropbyregion(L,DEM);
