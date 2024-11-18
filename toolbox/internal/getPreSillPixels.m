function PreSillPixel = getPreSillPixels(Z,FLATS,SILLS)
%GETPRESILLPIXELS Find pixels upstream of sill pixels
%
% Syntax
%
%     PreSillPixel = getPreSillPixels(Z,SILLS)
%
% Input arguments
%
%     Z      numeric matrix with elevations (filled DEM)
%     SILLS  Sill pixels (logical matrix)
%
% Output arguments
%
%     PreSillPixel Vector of linear indices of pixels 
%
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 06. September, 2024

siz = size(Z);

% establish the connectivity between sills and flats
[row,col] = find(SILLS);
IXsill    = sub2ind(siz,row,col);
rowadd    = [-1 -1 0 1 1  1  0 -1];
coladd    = [ 0  1 1 1 0 -1 -1 -1];
PreSillPixel = [];

for r = 1:8
    rowp = row + rowadd(r);
    colp = col + coladd(r);

    ValidRowColPair    = rowp>0 & colp>0 & rowp<=siz(1) & colp<=siz(2);
    IXPreSill = sub2ind(siz,rowp(ValidRowColPair),colp(ValidRowColPair));

    PreSillPixel = [PreSillPixel;...
        IXPreSill((Z(IXsill(ValidRowColPair)) == Z(IXPreSill)) & FLATS(IXPreSill))];   %#ok<AGROW>

end
