function [X,Y,HG] = hexgrid(DEM,d)

%HEXGRID creates an array of haxagonal points
%
% Syntax
%
%     [X,Y] = hexgrid(DEM)
%     [X,Y] = hexgrid(DEM,d)
%     [X,Y,HG] = ...
%
% Description
%
%     hexgrid creates an array of haxagonal points that are evenly
%     distributed across the entire DEM. The X and Y arrays are the center
%     coordinates of the grid cells. Currently only limited options.
%
% Input arguments
%
%     DEM   GRIDobj
%     d     factor (scalar) by which the grid cellsize is multiplied in the
%           x-direction to yield the desired point spacing. Default is 10.
%
% Output arguments
%
%     X,Y   x,y coordinate arrays
%     HG    geotable with hexagonal mappolyshapes with their centers at X
%           and Y (Note that a small value of d will slow down the function
%           substantially).
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     [X,Y] = hexgrid(DEM,100);
%     figure
%     imagesc(DEM)
%     hold on
%     scatter(X(:),Y(:),50,'ko','filled')
%     [VX,VY] = voronoi(X,Y);
%     plot(VX,VY,'w-','LineWidth',1)
%     axis(getextent(DEM))
%
% See also: meshgrid
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: 19. June 2025

arguments
    DEM  GRIDobj
    d {mustBePositive} = 10
end


[x,y] = getcoordinates(DEM);
cs = DEM.cellsize;

dx = cs*d;
dy = dx*sqrt(3)/2;

minx = min(x);
maxx = max(x);
miny = min(y);
maxy = max(y);

x = minx:dx:maxx;
y = miny:dy:maxy;

[X,Y] = meshgrid(x,y);

[nr,nc] = size(X);
sx = zeros(nr,1);
sx(1:2:end) = sx(1:2:end)+1;
dX = repmat(sx.*dx/2,[1,nc]);
X = X + dX;

%% Create hexgrid of map- or geopolyshapes
if nargout == 3
    % Voronoi
    % --- Input: Nx2 array of center points ---
    centers = [X(:) Y(:)];  % Replace with your actual center coordinates

    % --- Hexagon parameters ---
    r = d*cs / sqrt(3);  % Radius of hexagon (distance from center to any vertex)

    % Define hexagon shape around origin
    theta = (0:6) * pi/3 + pi/6;  % 6 vertices of regular hexagon
    x_hex = r * cos(theta);
    y_hex = r * sin(theta);

    % --- Create array of mappolyshape objects ---
    %numPoints = size(centers, 1);
    %polys(numPoints, 1) = mappolyshape();  % Preallocate

    x_shifted = x_hex + centers(:,1);
    y_shifted = y_hex + centers(:,2);

    x_shifted = num2cell(x_shifted,2);
    y_shifted = num2cell(y_shifted,2);

    if isGeographic(DEM)
        polys = cellfun(@(x,y)geopolyshape(x,y),y_shifted,x_shifted);
    else
        polys = cellfun(@(x,y)mappolyshape(x,y),x_shifted,y_shifted);
    end

    HG = table(polys,(1:numel(polys))',VariableNames= {'Shape','ID'});
    if isGeographic(DEM)
        HG.Shape.GeographicCRS = parseCRS(DEM);
    elseif isProjected(DEM)
        HG.Shape.ProjectedCRS = parseCRS(DEM);
    end
end
