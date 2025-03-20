function mustBeGRIDobjOrNalOrScalar(x,S)

%mustBeGRIDobjOrNalOrScalar Validate that value is GRIDobj, a valid nal or scalar
%
% Syntax
%
%     mustBeGRIDobjOrNalOrScalar(x,S)
%
% Description
%
%     The function validates that a value x is a GRIDobj or node-attribute
%     list that is aligned with a STREAMobj S or PPS object. x can also be
%     a scalar.
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 10. November, 2024

if isa(S,'PPS')
    S = S.S;
end

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

    TF = isnal(S,x) || isscalar(x);
    if ~TF  
        nnodes = numel(S.x);
        error('TopoToolbox:incorrectinput',...
            ['Input is not a valid node-attribute list of the supplied \n' ...
             'STREAMobj. To be valid, the input must be a column vector \n' ...
             'with ' num2str(nnodes) ' elements.']);

    end

end
