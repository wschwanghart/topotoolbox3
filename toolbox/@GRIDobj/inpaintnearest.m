function [A,D,IX] = inpaintnearest(A)

%INPAINTNEAREST Fill nan-pixels with values of nearest nan-pixels
%
% Syntax
%
%     B = inpaintnearest(A)
%     [B,D,L] = ...
%
% Description
%
%     inpaintnearest uses nearest neighbor interpolation to fill nan-pixels
%     in the GRIDobj A.
%
% Input arguments
%
%     A     GRIDobj
%
% Output arguments
%
%     B     GRIDobj with inpainted values
%     D     GRIDobj with distances
%     L     array with indices to non-nan pixels
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     % Set some values to nan
%     DEM = clip(DEM,DEM>1500 & DEM < 2000);
%     DEMi = inpaintnearest(DEM);
%     imageschs(DEMi)
%
%
% See also: GRIDobj, GRIDobj/inpaintnans, bwdist, bwdist_old
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 10. March, 2026

arguments
    A   GRIDobj
end

I     = isnan(A.Z);
% Due to a change in bwdist in 2025b which accepts no second output for
% arrays exceeding 2^24 elements, we revert to bwdist_old, if
% necessary.
if isMATLABReleaseOlderThan('R2025b') || numel(A.Z) < 2^24
    [D,L] = bwdist(~I,'euclidean');
else
    [D,L] = bwdist_old(~I,'euclidean');
end

A.Z = A.Z(L);

if nargout > 1
    D = GRIDobj(A,D*(A.cellsize));
end
if nargout > 2
    IX = L;
end