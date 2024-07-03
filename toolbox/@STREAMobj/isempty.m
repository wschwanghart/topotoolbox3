function tf = isempty(S)

%ISEMPTY Determine whether a STREAMobj is empty
%
% Syntax
%
%     tf = isempty(S)
%
% Description
%
%     isempty determines if S has at least two nodes and one edge. 
%
% Input arguments
%
%     S   STREAMobj
%
% Output arguments
%
%     tf  true or false
%
% See also: STREAMobj
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 2. July, 2024

tf = isempty(S.x);