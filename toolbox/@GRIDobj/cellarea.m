function CA = cellarea(DEM)

%CELLAREA Calculate cell areas of a GRIDobj in geographic coordinate system
%
% Syntax
%
%     CA = cellarea(DEM)
%
% Description
%
%     cellarea returns the area for each cell in the DEM with a geographic 
%     coordinate system (requires the mapping toolbox). DEM must have a
%     geographic coordinate system.
%
% Input arguments
%
%     DEM     GRIDobj
%
% Output arguments
%
%     CA      cell areas (GRIDobj) in m^2.
%
% See also: GRIDobj/isGeographic, GRIDobj/isProjected 
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 3. June, 2024

arguments
    DEM   GRIDobj
end

if isProjected(DEM)
    error("DEM must have a geographic coordinate system.")
end

if isGeographic(DEM)
    % There is a geocellreference or geopostingsreference
    unit = DEM.georef.GeographicCRS.Spheroid.LengthUnit;
    unit = validatestring(unit,{'meter','kilometer'});
    [~,ca] = areamat(true(DEM.size),DEM.georef);
    total_surface_area = 1;
    switch unit
        case 'kilometer'
            scale = 1e-6;
        otherwise
            scale = 1;
    end
else
    error('DEM must have a geographic coordinate system.')
end

ca = total_surface_area * ca * scale;

CA   = DEM;
CA.Z = repmat(ca,1,DEM.size(2));




