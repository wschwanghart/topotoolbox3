function [L,XLine,YLine,D2Line,T] = split2label(I,x,y)

%SPLIT2LABEL Split a GRIDobj into labels based on dividing line
%
% Syntax
%
%     L = split2label(I,x,y)
%     [L,XLine,YLine,D2Line,T] = ...
%
% Description
%
%     This function takes a GRIDobj I and returns a GRIDobj L where each
%     pixel is labelled according to its location relative to a line
%     defined by the coordinates of two points given by the 2x1 vectors x
%     and y. The function assumes an infinite line extending beyond the
%     line coordinates.
%
% Input arguments
%
%     I    GRIDobj
%     x,y  2x1 coordinate vectors that define a line
%
% Output arguments
%
%     L             labelled regions (-1 for pixels right to line, 1 for 
%                   pixels left to line, and 0 for pixels on line)
%     XLine,YLine   GRIDobjs of projected coordinates of each pixel on line
%     D2Line        Distance of each pixel to line. 
%     T             Relative distance along line. Pixels with values in the
%                   range of [0,1] are projected to locations on the line
%                   inbetween the coordinates in x and y.  
%
% Example 1
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     xl = [409959 376329]';
%     yl = [3800313 3792543]';
%     [L,~,~,D,T] = split2label(DEM,xl,yl);
%     subplot(1,2,1); imageschs(DEM,L); 
%     subplot(1,2,2); imageschs(DEM,D); 
%
% Example 2
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     tiledlayout(3,3,"TileSpacing","none","padding","tight")
%     for r = 1:9
%        [xl,yl] = randomsample(DEM,2);
%        L = split2label(DEM,xl,yl);
%        nexttile;
%        imagesc(L)
%     end
%
% See also: GRIDobj/dist2line, GRIDobj/dist2curve
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 29. January, 2026

arguments
    I GRIDobj
    x {mustBeNumeric}
    y {mustBeNumeric}
end

x = x(:);
y = y(:);

validateattributes(x,'numeric',{'numel',2},'split2label','x',2);
validateattributes(y,'numeric',{'numel',2},'split2label','y',2);

x1 = x(1);
x2 = x(2);
y1 = y(1);
y2 = y(2);

% Direction vector of line
dx = x2-x1;
dy = y2-y1;

% Raster coordinates
[xi,yi] = getcoordinates(I);

% Parameter t for projection
t = ((xi - x1)*dx + (yi - y1)*dy) / (dx^2 + dy^2);
 
% Signed cross product
s = (x2 - x1).*(yi - y1) - (y2 - y1).*(xi - x1);

% Prepare output
L = GRIDobj(I,reshape(sign(s),I.size));
if nargout > 1
    XLine = GRIDobj(I,x1+t.*dx);
    YLine = GRIDobj(I,y1+t.*dy);

    D2Line = GRIDobj(I,hypot(XLine.Z-xi,YLine.Z-yi));
    T  = GRIDobj(I,t);
end