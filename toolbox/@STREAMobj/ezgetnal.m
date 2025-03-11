function x = ezgetnal(S,x,cl)

%EZGETNAL Easy handling and retrieval of node-attribute lists
%
% Syntax
%
%     x = ezgetnal(S)
%     x = ezgetnal(S,val)
%     x = ezgetnal(S,val,class)
%
% Description
%
%     This is a small function that returns a node-attribute list of double
%     format, by default. Class can be any valid basic class (e.g.
%     'uint8').
%   
%     ezgetnal(S) returns a node-attribute list of zeros.
% 
%     ezgetnal(S,val) returns a node-attribute list of any input that the
%     function can digest. If val is a scalar, x is a node-attribute list 
%     where each element is val. If val is a GRIDobj, ezgetnal extracts the
%     values along the stream network. This is the same output as the
%     function getnal. If val is a node-attribute list, x equals val.
%     However, by default, x is converted to class double.
%
%     ezgetnal(S,val,class) lets you define the output class (e.g. 'double',
%     'single', 'logical'). The keyword 'same' will set the output class to
%     the input underlying class.
%
%   
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 17. July, 2024

arguments
    S    STREAMobj
    x    = 0
    cl   = 'double'
end

switch cl
    case 'same'
        cl = underlyingType(x);
end

if isa(x,'GRIDobj')
    % If x is a GRIDobj
    x = cast(getnal(S,x),cl);
elseif isnal(S,x)
    % If x is already a node-attribute list
    x = cast(x,cl);
elseif isscalar(x)
    % If x is a scalar, x will be implicitly expanded.
    x = cast(getnal(S) + x,cl);
else
    error('TopoToolbox:Input','Cannot handle input')
end
end

