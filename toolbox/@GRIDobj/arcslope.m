function SLP = arcslope(DEM,unit)

%ARCSLOPE Slope from a digital elevation model sensu ArcGIS
%
% Syntax
%
%     SLP = arcslope(DEM)
%     SLP = arcslope(DEM,unit)
%
% Description
%
%     ARCSLOPE returns the slope as calculated by ArcGIS (mean slope of
%     8 connected neighborhood). Default unit is as tangent, but you can 
%     specify alternative units identical to the gradient8 function.
%
% Input
%
%     DEM       digital elevation model (class: GRIDobj)
%     unit      'tan' --> tangent (default)
%               'rad' --> radian
%               'deg' --> degree
%               'sin' --> sine
%               'per' --> percent
% 
% Output
%
%     SLP         mean gradient (class: GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     SLP = arcslope(DEM,'deg');
%     imageschs(DEM,SLP);
%
% See also: GRIDobj/gradient8
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de) and
%         Adam M. Forte (aforte[a]asu.edu) 
% Date: 3. June, 2024

arguments
    DEM   GRIDobj
    unit  {mustBeTextScalar} = "tangent"
end

unit = validatestring(unit,{'tangent' 'degree' 'radian' 'percent' 'sine'},'arcslope','unit',2);

z=DEM.Z;

% Pad array to avoid edge effects
zp=padarray(z,[1 1],'replicate');

% Handle nans
I = isnan(zp);
in = any(I(:));
if in
    [~,L] = bwdist(~I);
    zp = zp(L);
end

% Define anon function to calculate the mean gradient of 8 connected neighborhood 
% by same algorithm as the ArcGIS slope function
kernel = [ -1 -2 -1; 0 0 0; 1 2 1]';
rr = sqrt((conv2(zp,kernel','valid')./(8*DEM.cellsize)).^2 + ...
                (conv2(zp,kernel,'valid')./(8*DEM.cellsize)).^2);

if in
    rr(I(2:end-1,2:end-1)) = nan;
end

% create a copy of the DEM instance
SLP = DEM;
SLP.name = 'slope (arcgis)';
SLP.zunit = unit;
SLP.Z=rr;

switch unit
    case 'tangent'
        % do nothing
    case 'degree'
        SLP.Z = atand(SLP.Z);
    case 'radian'
        SLP.Z = atan(SLP.Z);
    case 'sine'
        SLP.Z = sin(atan(SLP.Z));
    case 'percent'
        SLP.Z = SLP.Z*100;
end

