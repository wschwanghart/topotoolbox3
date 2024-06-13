function [S,locb] = trunk(S,A)

%TRUNK extract trunk stream (longest stream) 
%
% Syntax
%
%     S2 = trunk(S)
%     S2 = trunk(S,A)
%     S2 = trunk(S,a)
%     [S2,loc] = trunk(...)
%
% Description
%
%     trunk reduces a stream network to the longest streams in each stream
%     network tree (e.g. connected component). The algorithm identifies
%     the main trunk by sequently tracing the maximum downstream 
%     distance in upstream direction. 
%
%     If the second input argument is a flow accumulation grid A (or node
%     attribute list a) than the function will extract the main trunk by
%     tracing the maximum upstream area in upstream direction.
%
% Input 
%
%     S    stream network (STREAMobj)
%     A    GRIDobj with flow accumulation values (as returned by the
%          function flowacc). 
%     a    node-attribute list (nal) (e.g. getnal(S,flowacc(FD))).
% 
% Output
%
%     S2   stream network (STREAMobj) with only trunk streams in each
%          connecoted component
%     loc  linear index into node-attribute list of S
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     St = trunk(S);
%     plot(S)
%     hold on
%     plot(St)
%     legend('Stream network','Trunk rivers')
%
%
% See also: chiplot, FLOWobj/flowpathextract, STREAMobj/klargestconncomps,
%           STREAMobj/istrunk, STREAMobj/modify
%  
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 12. June, 2024


arguments
    S   STREAMobj
    A   {mustBeGRIDobjOrNalOrEmpty(A,S)} = []
end

narginchk(1,2);


if nargin == 2
    a = ezgetnal(S,A);
end

        
if nargout > 1
    % make copy
    Scopy = S;
end

% downstream distance
nrc = numel(S.x);

if nargin == 1
    dds = distance(S,'max_from_ch');
else
    dds = a;
end

D        = sparse(double(S.ix),double(S.ixc),dds(S.ix)+1,nrc,nrc);
OUTLET   = any(D,1)' & ~any(D,2);
[~,Imax] = max(D,[],1);
II       = false(nrc,1);
II(Imax) = true;

I = false(nrc,1);
I(OUTLET) = true;
for r = numel(S.ix):-1:1
    I(S.ix(r)) = I(S.ixc(r)) && II(S.ix(r));
end

L = I;
I = L(S.ixc) & L(S.ix);

S.ix  = S.ix(I);
S.ixc = S.ixc(I);

IX    = cumsum(L);
S.ix  = IX(S.ix);
S.ixc = IX(S.ixc);

S.x   = S.x(L);
S.y   = S.y(L);
S.IXgrid   = S.IXgrid(L);

if nargout > 1
    % Get indices to be able to index in node-attribute lists of Scopy
    [~,locb] = ismember(S.IXgrid,Scopy.IXgrid);
end
