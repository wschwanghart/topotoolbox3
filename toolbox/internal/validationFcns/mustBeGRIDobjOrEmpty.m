function mustBeGRIDobjOrEmpty(A)

%mustBeGRIDobjOrEmpty Validate that value is GRIDobj or an empty array
%
% Syntax
%
%     mustBeGRIDobjOrNal(x)
%
% Description
%
%     The function validates that a value x is a GRIDobj or an empty array
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 18. November, 2024

    if ~(isempty(A) || isa(A,'GRIDobj'))
        eid = 'Input:mustBeGRIDobjOrEmpty';
        msg = 'Input must be GRIDobj or empty.';
        error(eid,msg)
    end
end