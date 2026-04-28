function varargout = getoutline(DEM,options)

%GETOUTLINE Get outline of GRIDobj
%
% Syntax
%
%     GT = getoutline(DEM)
%     [x,y] = getoutline(DEM)
%     ... = getoutline(DEM,'pn',pv')
%
% Description
%
%     getoutline returns the outline of a GRIDobj. By default, getoutline
%     returns the coordinate vectors of the DEM edges. By setting shownans
%     = true, you can get the outline around the valid (non-nan) data in
%     the DEM. The function returns a geotable or nan-punctuated vectors of
%     the coordinates.
%
% Input arguments
%
%     DEM         GRIDobj
%     
%     Parameter name/value pairs
%     
%     'shownans'    {false} or true.
%
% Output arguments
%
%     GT     geotable that can be displayed using mapshow or geoplot. GT
%            can be exported with shapewrite.            
%     x,y    coordinate vectors that can be used to plot the extent.            
%
% Example 1
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     GT  = getoutline(DEM);
%     geoplot(GT)
%
% Example 2
%
%     DEM = readexample('taalvolcano');
%     I = identifyflats(DEM);
%     I.Z = bwareaopen(I.Z,20);
%     DEM.Z(I.Z) = nan;
%     GT = getoutline(DEM,'shownans',true);
%     geoplot(GT)
%     
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 20. April, 2026


arguments (Input)
    DEM  GRIDobj
    options.shownans = 0
end

% Check if there are nans
I = isnan(DEM);
% It's only required to account for nans if there are nans in the grid
nnan = options.shownans && any(I);

if ~nnan
    
    [x,y] = getcoordinates(DEM);

    csh   = DEM.cellsize/2;
    maxx = max(x)+csh;
    minx = min(x)-csh;    
    maxy = max(y)+csh;
    miny = min(y)-csh;
    
    x = [minx minx maxx maxx minx];
    y = [miny maxy maxy miny miny];

    
    if nargout == 1
    if isGeographic(DEM)
        GT = table(y,x);
        GT.Properties.VariableNames = {'Lat','Lon'};
        GT = table2geotable(GT,"CoordinateReferenceSystem",DEM.georef.GeographicCRS,...
            "GeometryType","polygon");
        GT = removevars(GT,["Lat" "Lon"]);
        
    elseif isProjected(DEM)
        GT = table(x,y);
        GT.Properties.VariableNames = {'X','Y'};
        GT = table2geotable(GT,"CoordinateReferenceSystem",DEM.georef.ProjectedCRS,...
            "GeometryType","polygon");
        GT = removevars(GT,["X" "Y"]);
    else
        GT = table(x,y);
        GT.Properties.VariableNames = {'X','Y'};
        GT = table2geotable(GT,"GeometryType","polygon");
        GT = removevars(GT,["X" "Y"]);

    end
    else
        x = x(:);
        y = y(:);
    end
    

else
    I = ~I;
    
    [GT,x,y] = GRIDobj2geotable(I,'conn',8,...
        'parallel',false);

end

if nargout == 1
    GT.ID = (1:size(GT,1))';
    varargout{1} = GT;

elseif nargout == 2
    % Two outputs, return x and y coordinates

    varargout{1} = x;
    varargout{2} = y;
end
end
