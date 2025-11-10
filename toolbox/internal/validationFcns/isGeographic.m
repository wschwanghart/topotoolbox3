function tf = isGeographic(A)
%ISGEOGRAPHIC Determines whether a TT object has a geographic coordinate system
%
% Syntax
%
%     tf = isgeographic(A)
%
tf = isprop(A.georef,"GeographicCRS");
if isempty(tf)
    tf = false;
    return
end
tf = tf && ~isempty(A.georef.GeographicCRS);
end