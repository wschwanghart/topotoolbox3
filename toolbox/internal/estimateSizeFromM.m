function siz = estimateSizeFromM(M)
%ESTIMATESIZEFROMM Estimate DEM size from sparse transfer matrix
%
% Syntax
%
%     siz = estimateSizeFromM(M)
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 06. September, 2024

[ic,icd] = find(M);
offsets  = unique(ic-icd);
siz = max(abs(offsets)-1);
siz(2) = size(M,1)./siz(1);

