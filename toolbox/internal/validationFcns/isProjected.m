function tf = isProjected(A)
%ISPROJECTED Determines whether GRIDobj has a projected coordinate system
%
% Syntax
%
%     tf = isProjected(A)
%

if isa(A,'PPS')
    A = A.S;
end

tf = isprop(A.georef,"ProjectedCRS");
if isempty(tf)
    tf = false;
    return
end
tf = tf && ~isempty(A.georef.ProjectedCRS);
end