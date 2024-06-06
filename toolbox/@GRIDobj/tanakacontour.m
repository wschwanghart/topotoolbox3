function tanakacontour(DEM,nlevels)

%TANAKACONTOUR Relief depiction using Tanaka contours
%
% Syntax
%
%     tanakacontour(DEM,nlevels)
%
% Description
%
%     Tanaka contours, also known as illuminated contours, are a
%     cartographic technique used to enhance the visual representation of
%     contour lines on maps. This method was developed by Japanese
%     geographer Kitiro Tanaka in the early 20th century. Tanaka contours
%     aim to simulate the effect of light and shadow to create a
%     three-dimensional appearance, making the terrain easier to interpret.
%
% Input arguments
%
%     DEM       GRIDobj
%     nlevels   number of contours (default = 10)
%
% See also: contourf, hillshade
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 6. June, 2024

arguments
    DEM  GRIDobj
    nlevels = 10
end

% calculate solar vector
azimuth = 315;
azid    = azimuth-90;
azisource = azid/180*pi;
[sx,sy,~] = sph2cart(azisource,0.1,1);

% get coordinate matrices and z values
[Z,X,Y] = GRIDobj2mat(DEM);

% calculate surface normals
[Nx,Ny] = surfnorm(Z);
N = hypot(Nx,Ny);
Nx = Nx./N;
Ny = Ny./N;

H = GRIDobj(DEM);
H = H-min(H);
H = H./max(H);
H.Z = Nx*sx + Ny*sy;


% Get coordinate matrices and Z
[~,h] = contourf(X,Y,Z,nlevels);
axis image
drawnow

% colormap
cmap = gray(255)*255;
cval = linspace(0,1,size(cmap,1))';

% the undocumented part
for r = 1:numel(h.EdgePrims)
    % the EdgePrims property contains vertex data
    x = h.EdgePrims(r).VertexData(1,:);
    y = h.EdgePrims(r).VertexData(2,:);
    % interpolate hillshade grayscale values to contour vertices
    c = interp(H,double(x(:)),double(y(:)));
    % interpolate to colormap
    clr = interp1(cval,cmap,c);
    % adjust color data to conform with EdgePrims ColorData property
    clr = clr';
    clr = uint8([clr;repmat(255,1,size(clr,2))]);
    % set all properties at once
    set(h.EdgePrims(r),'ColorBinding','interpolated','ColorData',clr,'LineWidth',1.2);
end
    



