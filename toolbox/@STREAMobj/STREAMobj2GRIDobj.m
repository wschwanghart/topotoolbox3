function G = STREAMobj2GRIDobj(S,nal)

%STREAMOBJ2GRIDOBJ Convert STREAMobj to GRIDobj
%
% Syntax 
%
%     G = STREAMobj2GRIDobj(S)
%     G = STREAMobj2GRIDobj(S,nal)
%
% Description
%     
%     STREAMobj2GRIDobj converts an instance of STREAMobj S to a new
%     instance of GRIDobj G with the same spatial reference and extent as
%     the instance from which S has been derived. If the second input is a
%     node-attribute list, STREAMobj2GRIDobj will write the node values to
%     the pixels in the output GRIDobj G.
%     
% Input arguments
%
%     S     STREAMobj
%     nal   node attribute list
%
% Output arguments
%
%     G     GRIDobj (contains logical raster with ones at stream locations
%           and zeros elsewhere. If a node attribute list was supplied, G
%           will have the same underlying class as nal. All non-stream
%           values are set to nan.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     g = gradient(S,DEM);
%     G = STREAMobj2GRIDobj(S,g);
%     imagesc(G)
%
% See also: STREAMobj, GRIDobj
% 
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 13. June, 2024

arguments
    S    STREAMobj
    nal  {mustBeGRIDobjOrNalOrEmpty(nal,S)} = [] 
end

if isempty(nal)
    nal = true(size(S.x));
    G = GRIDobj(S,'logical');
    G.Z(S.IXgrid) = nal;
else
    if isa(nal,'GRIDobj')
        cl = underlyingType(nal);
    else
        cl = class(nal);
    end
    nal = ezgetnal(S,nal,cl);
    G   = GRIDobj(S,'cl')*cast(nan,cl);
    G.Z(S.IXgrid) = nal;
end

G.name = 'stream grid';
