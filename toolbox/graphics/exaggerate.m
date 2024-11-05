function exaggerate(axes_handle,exagfactor)

%EXAGGERATE Elevation exaggeration in a 3D surface plot
%
% Syntax
%
%     exaggerate(axes_handle,exagfactor)
%
% Description
%
%     exaggerate is a simple wrapper for calling set(gca...). It controls
%     the data aspect ratio in a 3D plot and enables elevation
%     exaggeration.   
%
% Input
%
%     axes_handle   digital elevation model
%     exagfactor    exaggeration factor (default = 1)
%
% Example
%
%     load exampledem
%     for r = 1:4;
%     subplot(2,2,r);
%     surf(X,Y,dem); exaggerate(gca,r);
%     title(['exaggeration factor = ' num2str(r)]);
%     end
% 
% See also: GRIDobj/SURF, exaggerate, GRIDobj/imageschs
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 4. November 2024

arguments
    axes_handle {mustBeA(axes_handle,'matlab.graphics.axis.Axes')} = gca
    exagfactor (1,1) {mustBeNumeric} = 1
end

axis(axes_handle,'image');
set(axes_handle,'DataAspectRatio',[1 1 1/exagfactor]);
