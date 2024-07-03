function S = klargestconncomps(S,k)

%KLARGESTCONNCOMPS Keep k largest connected components in an instance of STREAMobj
%
% Syntax
%
%     S2 = klargestconncomps(S,k)
%
% Description
%
%     klargestconncomps extracts the k largest strongly connected 
%     components from the stream network in an instance of STREAMobj.
%     Largest refers to the number of edges (links between individual
%     pixels) in the network.
%
% Input arguments
%
%     S        instance of STREAMobj
%     k        number of largest components (default k=1)
%
% Output arguments
%
%     S2       new instance of STREAMobj that contains only the k largest
%              components of S
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     S  = STREAMobj(FD,A>1000);
%     S2 = klargestconncomps(S,2);
%     plot(S)
%     hold on
%     plot(S2)
%     legend('Stream network','2 largest conn. components')
%
% See also: STREAMobj/modify, STREAMobj/conncomps
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 12. June, 2024


arguments
    S    STREAMobj
    k    (1,1) {mustBePositive,mustBeInteger} = 1
end

% create sparse adjacency matrix
nrc = numel(S.x);
M = sparse(double(S.ix),double(S.ixc),true,nrc,nrc);
M = M | M' | speye(nrc);
% use dmperm to find connected components
% see http://www.mathworks.de/matlabcentral/fileexchange/21366-findcomponents
% by Tim Davis
[~,p,~,r] = dmperm(M);
nc = length(r) - 1;

% label matrix
L = zeros(nrc,1);
% sort to find the largest conncompoments
[~,dd] = sort(diff(r),'descend');
if k>nc
    warning('TopoToolbox:STREAMobj',...
            ['There are only ' num2str(nc) ' connected components in the stream network']);
end

k = min(nc,k);
for tt = dd(1:k) 
    L(p(r(tt):r(tt+1)-1)) = tt;
end

S = subgraph(S,L>0);

