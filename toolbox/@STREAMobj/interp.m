function zi = interp(S,d,z,method,extrapolation)

%INTERP Interpolate data on STREAMobj (single river only)
%
% Syntax
%
%     zi = interp(S,d,z)
%     zi = interp(S,d,z,'method','extrapolation')
%
% Description
%
%     
% Input arguments
%
%     S      STREAMobj with one channelhead
%     d      distance values
%     z      variable values (same size as d)
%     
%     'method'  see interp1 for details
%     'extrapolation'  see interp1 for details
%
% Output arguments
%
%     zi     node-attribute list with interpolated values
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S   = STREAMobj(FD,'minarea',1000);
%     S   = trunk(klargestconncomps(S));
%     zi = interp(S,[10000 12000 20000 39000 40000],[5 10 20 6 10],...
%                 'spline',10);
%     stackedplotdz(S,{DEM zi})
%     
% See also: STREAMobj, STREAMobj/trunk, demo_modifystreamnet
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 11. June, 2026

arguments
    S {mustBeSingleRiver}
    d 
    z {mustHaveEqualSize(z,d)}
    method = 'linear'
    extrapolation = 'extrap'
end


[d,ix] = sort(d,'ascend');
z = z(ix);

zi = interp1(d,z,S.distance,method,extrapolation);

end

function mustBeSingleRiver(S)
if numel(streampoi(S,'channelheads','ix')) ~= 1
    error('STREAMobj must contain a single stream.')
end
end

function mustHaveEqualSize(x,y)
if numel(x) ~= numel(y)
    error('d and z must have the same number of elements.')
end
end
