function tf = isEqualGeoreference(A,B)
%ISEQUALGEOREFERENCE Determines whether two TT objects have the same georeferencing
%
% Syntax
%
%     tf = isEqualGeoreference(A,B)
%
% If A and B have no referencing objects, their reference is
% assumed to be the same.

% tf = false;


% Both GRIDobjs do not have referencing objects
if isempty(isGeographic(A)) && isempty(isGeographic(B))
    tf = true;
    return
end

% One of both GRIDobj does not have a referencing object
if xor(isempty(isGeographic(A)),isempty(isGeographic(B)))
    tf = false;
    return
end

% Both have different types of referencing objects
if xor(isGeographic(A),isGeographic(B))
    tf = false;
end

% They have geographic referencing objects
if isGeographic(A)
    tf = isequal(A.georef.GeographicCRS,...
        B.georef.GeographicCRS);
    return
end

% They have projected referencing objects
if isProjected(A)
    tf = isequal(A.georef.ProjectedCRS,...
        B.georef.ProjectedCRS);
    return
end
end
