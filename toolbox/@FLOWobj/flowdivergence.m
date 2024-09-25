function D = flowdivergence(FD)

%FLOWDIRVERGENCE Calculate the number of downstream neighbors
%
% Syntax
%
%     D = flowdivergence(FD)
%
% Description
%
%     flowdivergence determines the number of downstream neighbor of each
%     cell. 
%
% Input arguments
%
%     FD   FLOWobj
%
% Output arguments
%
%     D    GRIDobj 
%
%
% See also: FLOWobj, FLOWobj/ismulti
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 25. September, 2024

arguments
    FD    FLOWobj
end

D = accumarray(FD.ix,FD.ixc,[prod(FD.size) 1],@(x) uint8(numel(x)),zeros(1,'uint8'));
D = reshape(D,FD.size);
D = GRIDobj(FD,D);
