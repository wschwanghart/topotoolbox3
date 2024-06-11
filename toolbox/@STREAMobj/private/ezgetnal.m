function x = ezgetnal(S,x,cl)

%EZGETNAL Easy handling and retrieval of node-attribute lists
%
% Syntax
%
%     x = ezgetnal(S,x)
%     x = ezgetnal(S,x,class)
%
% Description
%
%     This is a small function that returns a node-attribute list of double
%     format, by default. Class can be any valid basic class (e.g.
%     'uint8').
%   
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 6. June, 2024

arguments
    S    STREAMobj
    x    
    cl   = 'double'
end

if isa(x,'GRIDobj')
    % If x is a GRIDobj
    x = cast(getnal(S,x),cl);
elseif isnal(S,x)
    % If x is already a node-attribute list
    x = cast(x,cl);
elseif isnumeric(x) && isscalar(x)
    % If x is a scalar, x will be implicitly expanded.
    x = cast(getnal(S) + x,cl);
else
    error('TopoToolbox:Input','Cannot handle input')
end
end

