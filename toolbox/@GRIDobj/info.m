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
% Output arguments
%
%     -
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     info(DEM)
%
% See also: DISP
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 17. October, 2024

arguments
    DEM  GRIDobj
end

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
disp(['  data type:             ' underlyingType(DEM)])
disp(['  number of rows:        ' num2str(DEM.size(1))])
disp(['  number of columns:     ' num2str(DEM.size(2))])
disp(['  cellsize:              ' num2str(DEM.cellsize)])
disp(['  extent in x-direction: ' num2str(ul(1)) ' -- ' num2str(lr(1))])
disp(['  extent in y-direction: ' num2str(ul(2)) ' -- ' num2str(lr(2))])
disp(['  maximum z-value:       ' num2str(max(DEM))])
disp(['  minimum z-value:       ' num2str(min(DEM))])
disp(['  z-unit:                ' DEM.zunit])

if isGeographic(DEM)
disp(['  coordinate system:     ' ' Geographic coordinate system'])
disp(['                          '  char(DEM.georef.GeographicCRS.Name)])
elseif isProjected(DEM)
disp(['  coordinate system:     ' ' Projected coordinate system'])
disp(['                          '  char(DEM.georef.ProjectedCRS.Name)])
else
disp(['  coordinate system:     ' ' undefined (.georef empty)'])
end

    
disp(' ')


