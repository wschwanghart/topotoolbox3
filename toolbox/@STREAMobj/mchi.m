function mc = mchi(S,DEM,chi)

%MCHI Gradient of stream profile in chi space (M_chi = Ksn)
%
% Syntax
%
%     mc = mchi(S,DEM,chi)
%
% Description
%
%     mchi calculates the gradient of elevation in a DEM where the
%     horizontal coordinate is chi as obtained from the function
%     chitransform. mchi provides a metric that can be used to compare
%     channel gradients with different drainage areas and is the same as
%     ksn.
%
% Input arguments
%
%     S        STREAMobj
%     DEM      digital elevation model (GRIDobj)
%     z        node attribute list with elevation values
%     chi      node attribute list with chi values (see chitransform)
%
% Output arguments
%
%     mc       node attribute list with M_chi values
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1e5,'unit','map');
%     A = flowacc(FD);
%     c = chitransform(S,A,'mn',0.39);
%     z = imposemin(S,DEM);
%     mc = mchi(S,z,c);
%
%     % plot results
%     imageschs(DEM,DEM,'colormap',[1 1 1],'colorbar',false)
%     hold on
%     plotc(S,smooth(S,mc,'k',100));
%     h = colorbar;
%     h.Label.String = 'M_\chi';
%
% See also: STREAMobj/gradient, STREAMobj/chitransform, STREAMobj/ksn
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 2. July, 2024

arguments
    S  STREAMobj
    DEM {mustBeGRIDobjOrNal(DEM,S)}
    chi {mustBeNumeric}
end

% get node attribute list with elevation values
z = ezgetnal(S,DEM);

mc = zeros(size(z));
mc(S.ix) = (z(S.ix)-z(S.ixc))./(chi(S.ix)-chi(S.ixc));