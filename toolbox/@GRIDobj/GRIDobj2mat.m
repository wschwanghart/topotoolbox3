function [Z,x,y] = GRIDobj2mat(DEM)

%GRIDobj2mat Convert GRIDobj to matrix and coordinate vectors
%
% Syntax
%
%     [Z,X,Y] = GRIDobj2mat(DEM)
%
% Description
%
%     convert GRIDobj to matrix and coordinate vectors 
%
% Input
%
%     DEM       instance of GRIDobj class
% 
% Output
%
%     Z         matrix 
%     x,y       coordinate vectors
%                  
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     [Z,x,y] = GRIDobj2mat(DEM);
%     plot(x,Z(20,:))
%     xlabel('x-coordinate [m]')
%     ylabel('Elevation [m]')
%
% See also: GRIDobj 
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 3. June, 2024

Z = DEM.Z;
[x,y] = wf2XY(DEM.wf,DEM.size);
