function mustBeGRIDobjOrNal(x,S)

%mustBeGRIDobjOrNal Validate that value is GRIDobj or a valid nal
%
% Syntax
%
%     mustBeGRIDobjOrNal(x,S)
%
% Description
%
%     The function validates that a value x is a GRIDobj or node-attribute
%     list that is aligned with a STREAMobj S.
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 11. June, 2024

if isa(x,'GRIDobj')
    TF = isequal(S.size,x.size) && isequal(S.wf,x.wf);
    if ~TF     
        error('TopoToolbox:incorrectinput',...
            ['STREAMobj and GRIDobj do not align each other. Make sure that \n' ...
            'both instances have the same spatial reference. Both variables \n' ...
            'are deemed to have the same reference if their properties ''size'' \n' ...
            'and ''wf'' are both equal.']);

    end
else
    TF = isnal(S,x);
    if ~TF     
        nnodes = numel(S.x);
        error('TopoToolbox:incorrectinput',...
            ['Input is not a valid node-attribute list of the supplied \n' ...
             'STREAMobj. To be valid, the input must be a column vector \n' ...
             'with ' num2str(nnodes) ' elements.']);

    end

end
