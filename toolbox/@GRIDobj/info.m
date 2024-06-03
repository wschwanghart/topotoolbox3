function info(DEM)

%INFO Detailed information on GRIDobj instance
%
% Syntax
%
%     info(DEM)
%
% Description
%
%     info displays detailed information about an instance of an GRIDobj in
%     the command window.
%
% Input arguments
%
%     DEM    instance of GRIDobj
%
% See also: DISP
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 8. January, 2013

wf = DEM.wf;
s = inputname(1);
[ul] = wf*[1 1 1]';
[lr] = wf*[DEM.size 1]';

disp(' ')
if isempty(s)
    disp('TopoToolbox GRIDobj')
else
    disp(['TopoToolbox GRIDobj ' s ' (' ...
        '<a href = "matlab:openvar ' s '">show</a>/<a href = "matlab:imagesc(' s ')">plot</a>)']);
end
disp(['  name:                  ' DEM.name])
disp(['  data type:             ' class(DEM.Z)])
disp(['  number of rows:        ' num2str(DEM.size(1))])
disp(['  number of columns:     ' num2str(DEM.size(2))])
disp(['  cellsize:              ' num2str(DEM.cellsize)])
disp(['  extent in x-direction: ' num2str(ul(1)) ' -- ' num2str(lr(1))])
disp(['  extent in y-direction: ' num2str(ul(2)) ' -- ' num2str(lr(2))])
disp(['  maximum z-value:       ' num2str(max(DEM.Z(:)))])
disp(['  minimum z-value:       ' num2str(min(DEM.Z(:)))])
disp(['  z-unit:                ' DEM.zunit])

tf = isGeographic(DEM);
if isempty(tf)
    
disp(['  coordinate system:     ' ' undefined (.georef empty)'])
else
if tf
disp(['  coordinate system:     ' ' Geographic coordinate system'])
disp(['                          '  DEM.georef.GeographicCRS.Name])
else
disp(['  coordinate system:     ' ' Projected coordinate system'])
disp(['                          '  char(DEM.georef.ProjectedCRS.Name)])
end
end
    
disp(' ')


