function Lgrid = shufflelabel(Lgrid,r)

%SHUFFLELABEL Shufflelabel randomly relabels a label matrix
%
% Syntax
%
%     L = shufflelabel(L)
%     L = shufflelabel(L,reset)
%
% Description
%
%     shufflelabel randomly changes the order of labels in the label
%     GRIDobj L. Zeros and nans are ignored. The function is helpful when
%     plotting drainage basins.
%
%     When called with two input arguments, reset must be either true or 
%     false. If true, label values in L are reset to range from
%     one to numel(unique(L(:)). Note that this only works for 
%     numeric arrays.
%
% Input arguments
%
%     L      grid with ordinal or categorical data (class: GRIDobj)
%     reset  false (default) or true
%
% Output arguments
%
%     Ls     grid with values randomly shuffled
%
% Example
% 
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     L   = reclassify(DEM,'equalquantiles',10);
%     L   = shufflelabel(L);
%     imagesc(L)
%
%
% See also: RANDPERM, BWLABEL, LABELMATRIX
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 18. June, 2024


arguments
    Lgrid  GRIDobj
    r  (1,1) = false
end

% get values in GRIDobj
L = Lgrid.Z;

% is L a numeric array? 
inum = isnumeric(L);

% applies only if L is a numeric array.
validateattributes(r,{'numeric','logical'},{'scalar'})
r = r && inum;

% size of L
siz = size(L);

% force column vector
L   = L(:);

% exclude zeros and nans from shuffling (only when L is 
% a numeric array)
if inum
    I   = ~(L==0 | isnan(L));
else
    I   = ':';
end

% find unique elements in label vector
[uniqueL,~,ix] = unique(L(I));

% shuffle labels
if r
    uniqueLS = randperm(numel(uniqueL));
else
    uniqueLS = uniqueL(randperm(numel(uniqueL)));
end

% and map labels back into L
L(I) = uniqueLS(ix);

% finally reshape back to original size
L = reshape(L,siz);

% write to GRIDobj
Lgrid.Z = L;

