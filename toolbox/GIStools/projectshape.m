function GT = projectshape(GT,TARGET,projsource)

% PROJECTSHAPE Changes the projection of a geographic data structure
%
% Syntax
%
%     GT2 = projectshape(MS,TARGET)
%     GT2 = projectshape(MS,TARGET,projsource)
%     GT2 = projectshape(GT,TARGET)
%
% Description
%     
%     projectshape(MS,GRID) transforms the lat,lon of the features in the
%     geographic data structure MS (as obtained from importing a shapefile 
%     using shaperead) to the x,y units of the GRIDobj GRID.
%
%     projectshape(MS,GRID,mstruct) uses the map projection structure
%     mstruct to first convert the x,y units of the features in the
%     geographic data structure MS (as obtained from importing a shapefile 
%     using shaperead) before converting to the x,y units of the GRIDobj 
%     GRID.
% 
%     Note that projectshape requires the Mapping Toolbox.
%
% Input arguments
%
%     MS            geographic data structure
%     GT            geotable
%     GRID          GRIDobj, e.g., a DEM
%     projsource    projection (e.g. projcrs or geocrs) that is readable by
%                   the function parseCRS
%     
% Output arguments
%
%     GT2     geotable
%
% See also: parseCRS, geotable2mapstruct, mapstruct2geotable
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de) and
%         Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 16. June, 2024

arguments
    GT
    TARGET
    projsource = 4326
end

if isstruct(GT)
    GT = mapstruct2geotable(GT,'CoordinateReferenceSystem',projsource);
end

% Check input geometry type
geomType = GT.Shape.Geometry;

% Get source CRS
[CRSsource,issourceproj] = parseCRS(GT);
% Get target CRS
[CRStarget,istargetproj] = parseCRS(TARGET);

% Are both equal?
if isequal(CRSsource,CRStarget)
    return
end

% Convert geotable to table
if ~issourceproj
    T = geotable2table(GT, ["Latitude" "Longitude"]);
    if iscell(T.Latitude) % Input is polygon or line
        latc = T.Latitude;
        lonc = T.Longitude;
        ispoint = false;
    else % Input is point
        latc = num2cell(T.Latitude);
        lonc = num2cell(T.Longitude);
        ispoint = true;
    end
else
    T = geotable2table(GT, ["x" "y"]);
    if iscell(T.x) % Input is polygon or line
        xc = T.x;
        yc = T.y;
        ispoint = false;
    else % Input is point
        xc = num2cell(T.x);
        yc = num2cell(T.y);
        ispoint = true;
    end
end

% Now, all coordinates are available as cell arrays

if istargetproj && issourceproj
    [latc,lonc] = cellfun(@(x,y) projinv(CRSsource,x,y),xc,yc,...
        "UniformOutput",false);
    [xc,yc] = cellfun(@(lat,lon) projfwd(CRStarget,lat,lon),latc,lonc,...
        "UniformOutput",false);
elseif istargetproj && ~issourceproj
    [xc,yc] = cellfun(@(lat,lon) projfwd(CRStarget,lat,lon),latc,lonc,...
        "UniformOutput",false);
elseif issourceproj && ~istargetproj
    [latc,lonc] = cellfun(@(x,y) projinv(CRSsource,x,y),xc,yc,...
        "UniformOutput",false);
else
    error("Geocrs to geocrs is not supported.")
end

% Now create a new geotable
% If input geometry is point, then create numeric vectors
if ispoint && ~istargetproj
    latc = vertcat(latc{:});
    lonc = vertcat(lonc{:});
elseif ispoint && istargetproj
    xc = vertcat(xc{:});
    yc = vertcat(yc{:});
end
    
if ~istargetproj
    T.Latitude = latc;
    T.Longitude = lonc;

    T = table2geotable(T,'geographic',["Latitude" "Longitude"],...
        'CoordinateReferenceSystem',CRStarget,'GeometryType',geomType);
    
else
    T.x = xc;
    T.y = yc;

    T = table2geotable(T,'planar',["x" "y"],...
        'CoordinateReferenceSystem',CRStarget,'GeometryType',geomType);
end

GT.Shape = T.Shape;
