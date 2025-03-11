function z = cummaxupstream(S,z,cl)

%CUMMAXUPSTREAM Cumulative maximum in upstream direction
%
% Syntax
%
%     zm = cummaxupstream(S,DEM)
%     zm = cummaxupstream(S,z)
%     zm = cummaxupstream(S,z,class)
%
% Description
%
%     CUMMAXUPSTREAM calculates the cumulative maximum in upstream
%     direction. A value in zm will have the maximum value of its
%     downstream pixels in DEM or z. 
%
% Input arguments
%
%     S       STREAMobj
%     DEM     GRIDobj
%     z       node-attribute list
%     class   Output class. Default is 'same' (see ezgetnal). 
%
% Output arguments
%
%     zm    node-attribute list
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     S = trunk(S);
%     zm = cummaxupstream(S,DEM);
%     plotdz(S,DEM)
%     hold on
%     plotdz(S,zm)
%     hold off
%
% See also: STREAMobj/imposemin, cummax
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 7. March, 2025

arguments
    S  STREAMobj
    z {mustBeGRIDobjOrNal(z,S)}
    cl = 'same'
end

z = ezgetnal(S,z,cl);

ix = S.ix;
ixc = S.ixc;
for r = numel(ix):-1:1
    z(ix(r)) = max(z(ixc(r)),z(ix(r)));
end

