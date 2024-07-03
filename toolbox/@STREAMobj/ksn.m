function k = ksn(S,DEM,A,theta,K)

%KSN normalized steepness index
%
% Syntax
%
%     k = ksn(S,DEM,A)
%     k = ksn(S,z,a)
%     k = ksn(...,theta)
%     k = ksn(...,theta,K)
%
% Description
%
%     KSN returns the normalized steepness index using a default concavity
%     index of 0.45. The Ksn index, or normalized steepness index, is a
%     metric used in topographic and tectonic analyses to quantify the
%     steepness of river channels relative to a normalized reference
%     concavity. Steeper and possibly more erosive segments of the river
%     have higher Ksn values. The index helps identify areas of
%     differential uplift, varying lithologies, or tectonic activity by
%     comparing the steepness of river profiles across different regions.
%
% Input arguments
%
%     S      STREAMobj
%     DEM    digital elevation model (GRIDobj)
%     A      flow accumulation as returned by flowacc (GRIDobj). Note that
%            flowacc returns the number of pixels. The function ksn
%            calculates area in m^2 from these values internally.
%     z      node attribute list of elevation values
%     a      node attribute list of flow accumulation values
%     theta  concavity (default 0.45)
%     K      smoothing factor K (by default 0 = no smoothing). See function
%            STREAMobj/smooth for details. If K > 0, then the function will
%            smooth the river profile before river gradient and ksn are
%            calculated. To have more control on profile smoothing, see the
%            functions STREAMobj/crs (or crsapp) or STREAMobj/smooth.
%
% Output arguments
%
%     k      normalized steepness index
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     DEM = imposemin(S,DEM);
%     A = flowacc(FD);
%     k = ksn(S,DEM,A,0.45,100);
%     subplot(2,1,1);
%     imageschs(DEM,DEM,'colormap',[.9 .9 .9],'colorbar',false);
%     hold on
%     plotc(S,k)
%     
%     subplot(2,1,2);
%     plotdz(S,DEM,'color',k);
%
% See also: STREAMobj/crs, STREAMobj/smooth, FLOWobj/flowacc
% 
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 2. June, 2024

arguments
    S   STREAMobj
    DEM {mustBeGRIDobjOrNal(DEM,S)}
    A   {mustBeGRIDobjOrNal(A,S)}
    theta (1,1) {mustBeNumeric,mustBePositive} = 0.45
    K   (1,1) {mustBeNumeric,mustBeNonnegative} = 0
end

z = ezgetnal(S,DEM);
a = ezgetnal(S,A);

% minima imposition to avoid negative gradients
z = imposemin(S,z,0.00001);

% Smoothing, if required
if K ~= 0
    z = smooth(S,z,'K',K);
end

% calculate gradient
g = gradient(S,z);
% upslope area
% if ~isgeographic(S)
a = a.*S.cellsize.^2;
% end

k = g./(a.^(-theta));




