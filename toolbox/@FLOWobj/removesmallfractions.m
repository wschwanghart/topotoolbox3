function FD = removesmallfractions(FD,minfrac)

%removesmallfractions Remove links with small fractions in FLOWobj
%
% Syntax
%
%     FD = removesmallfractions(FD,minfrac)
%
% Description
%     
%     The function removes links in a FLOWobj where fractions indicate only
%     a small transfer amount to neighboring pixels. The operation makes
%     flow less dispersive.
%
% Input arguments
%
%     FD        FLOWobj
%     minfrac   positive scalar indicating the minimum fraction to remain
%               in the FLOWobj (must be less than 1/8).
%
% Output arguments
%
%     FD        FLOWobj
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEM = filter(DEM);
%     C = curvature(DEM,'planc');
%     imageschs(DEM,C,'percentclip',0.1)
%
% See also: FLOWobj/multi_normalize, FLOWobj/randomize
%        
% Author:  Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 2. January, 2025

arguments
    FD  FLOWobj
    minfrac (1,1) {mustBeLessThan(minfrac,0.1250),mustBePositive} = 0.1
end

switch FD.type
    case 'single'
        return
end

I = FD.fraction <= minfrac;
FD.ix(I) = [];
FD.ixc(I) = [];
FD.fraction(I) = [];
FD = multi_normalize(FD);

