function M = FLOWobj2M(FD)

%FLOWOBJ2M Convert instance of FLOWobj to flow direction matrix 
%
% Syntax
%
%     M = FLOWobj2M(FD);
%
% Description
%
%     FLOWobj2M converts an instance of FLOWobj to the flow direction
%     matrix M as used in TopoToolbox 1.
%
% Input arguments 
%
%     FD     FLOWobj
%
% Output arguments
%
%     M      sparse transfer matrix
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     M = FLOWobj2M(FD);
%     [x,y] = getcoordinates(DEM);
%     [x,y] = meshgrid(x,y);
%     gplot(M,[x(:) y(:)]) 
%
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 31. August, 2024

arguments 
    FD   FLOWobj
end

nrc = prod(FD.size);
switch lower(FD.type)
    case {'multi','dinf'}
        M = sparse(double(FD.ix),double(FD.ixc),FD.fraction,nrc,nrc);
    case 'single'
        M = sparse(double(FD.ix),double(FD.ixc),1,nrc,nrc);
end
end