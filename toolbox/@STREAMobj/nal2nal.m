function nal2 = nal2nal(S2,S1,nal1,fillval)

%NAL2NAL Map one node-attribute list to another 
%
% Syntax
%
%     nalB = nal2nal(SB,SA,nalA)
%     nalB = nal2nal(SB,SA,nalA,fillval)
%     nalB = nal2nal(SB,SA,nalA,nalB)
%
% Description
%
%     nal2nal maps the node-attribute list (nal) nalA of the stream network 
%     SA to another stream network SB. SA must be a subgraph of SB, i.e.,
%     SA is formed from a subset of the vertices (nodes) in SB.
%
% Input arguments
%
%     SB       STREAMobj
%     SA       STREAMobj that is a subgraph of SB
%     nalA     node-attribute list of SA
%     fillval  value to be assigned to nodes in nalB that are not members
%              of the network SA. By default, the value is NaN.
%     nalB     nal of SB. nal2nal will overwrite values at vertex-locations
%              in SA.
%     
% Output arguments
%
%     nalB     nal of SB
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S   = STREAMobj(FD,'minarea',1000);
%     S   = klargestconncomps(S);
% 
%     A   = flowacc(FD);
%     c   = chitransform(S,A); 
%     S2  = modify(S,'streamorder',1);
%     c2  = nal2nal(S2,S,c,0);
%     plot(S,'k')
%     hold on
%     plotc(S2,c2); 
%     h = colorbar;
%     h.Label.String ='\chi'; 
%
%
% See also: STREAMobj2GRIDobj
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]uni-potsdam.de)
% Date: 17. June, 2024

arguments
    S2    STREAMobj
    S1    STREAMobj
    nal1  {mustBeGRIDobjOrNal(nal1,S1)}
    fillval = nan
end

% Handle input nal
if isa(nal1,'GRIDobj')
    nal1 = ezgetnal(S,nal1,underlyingType(nal1));
end

% Handle output nal
if ~islogical(nal1)
    nal2 = ezgetnal(S2,fillval,class(nal1));
else
    if isscalar(fillval)
        fillval = false;
        nal2 = ezgetnal(S2,fillval,'logical');
    else
        nal2 = ezgetnal(S2,fillval,class(nal1));
    end 
end
    
[I,locb] = ismember(S1.IXgrid,S2.IXgrid);
nal2(locb(I)) = nal1(I);