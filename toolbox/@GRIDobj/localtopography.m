function OUT = localtopography(DEM,radius,options)

%LOCALTOPOGRAPHY Local topography
%
% Syntax
%
%     H = localtopography(DEM)
%     H = localtopography(DEM,radius)
%     H = localtopography(DEM,radius,pn,pv,...)
%
% Description
%
%     localtopography quantifies local relief, e.g. the elevation range
%     within a specific radius. localtopography may take a while to
%     evaluate large DEMs with large kernels. You may speed up calculations
%     by adjusting the parameter 'N' for calculating the max, min or range
%     filter. 
%
%     localtopography uses symmetric boundary padding. NaNs
%     are inpainted by nearest neighbor interpolation prior to the
%     calculation (affects only 'mean', 'median', 'prctile', 'std'). For
%     'max', 'min' and 'range', localtopography adopts the padding behavior
%     of the function imdilate and imerode.
%
% Input arguments
%
%     DEM     Digital elevation model (GRIDobj)
%     radius  radius of the moving window filter in map units. The default
%             value is 5000 (m).
%
% Parameter Name/Value (pn,pv) pairs
%
%     'type'   'range' (default), 'max', 'min', 'mean', 'median',
%              'prctile', 'std' (standard deviation)
%     'prc'    scalar between 0 and 100 [%]. Only applicable if 'prctile'
%              is chosen as 'type'.
%     'N'      enable speed improvement for 'max', 'min' and 'range'.
%              N must be 0, 4, 6,or 8 (default). When N is greater than 0,
%              the disk-shaped structuring element is approximated by a
%              sequence of N periodic-line structuring elements. When N
%              equals 0, no approximation is used, and the structuring
%              element members consist of all pixels whose centers are no
%              greater than R away from the origin. Choosing N entails a
%              trade-off between how well the shape of the disc-shaped
%              moving window is approximated and computational efficiency.
%              If you choose a large radius, then N=0 can take a while to 
%              evaluate.
%
% Output arguments
%
%     H      local topography grid (GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     H = localtopography(DEM,500);
%     imageschs(DEM,H)
%
% See also: IMDILATE, IMERODE
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 28. January, 2024

arguments
    DEM  GRIDobj
    radius (1,1) = 5000
    options.type {mustBeMember(options.type, ...
        {'range','max','min','mean','median','prctile','std'})} = 'range'
    options.N (1,1) {mustBeMember(options.N,[0 4 6 8])} = 8
    options.prc (1,1) {mustBeInRange(options.prc,0, 100,'exclusive')} = 90
end

dem = DEM.Z;
cs  = DEM.cellsize;

% any nans
INAN = isnan(dem);
flaginan = any(INAN(:));

% structuring element
radiuspx = ceil(radius/cs);
SE = strel('disk',radiuspx,options.N);

switch options.type
    case 'max'
        % Maximum filter
        if flaginan
            dem(INAN) = -inf;
        end
        H = imdilate(dem,SE);
    case 'min'
        % Minimum filter
        if flaginan
            dem(INAN) = inf;
        end
        H = imerode(dem,SE);
    case 'range'
        if flaginan
            dem(INAN) = -inf;
        end
        H1 = imdilate(dem,SE);
        if flaginan
            dem(INAN) = inf;
        end
        H2 = imerode(dem,SE);
        H  = H1-H2;
    case {'mean','average'}
        if flaginan
            DEM = inpaintnearest(DEM);
            dem = DEM.Z;
        end            
        H   = fspecial('disk',radiuspx);
        H   = imfilter(dem,H,'symmetric','same','conv');
    case 'median'
        if flaginan
            DEM = inpaintnearest(DEM);
            dem = DEM.Z;
        end
        H   = getnhood(SE);
        n   = round(sum(H(:))/2);
        H   = ordfilt2(dem,n,H,'symmetric');
        
    case 'prctile'
        if flaginan
            DEM = inpaintnearest(DEM);
            dem = DEM.Z;
        end
        H   = getnhood(SE);
        n   = round(sum(H(:))*options.prc/100);
        H   = ordfilt2(dem,n,H,'symmetric');
    
    case 'std'
        if flaginan
            DEM = inpaintnearest(DEM);
            dem = DEM.Z;
        end
        H   = getnhood(SE);
        H   = stdfilt(dem,H);        
end
            

% Handle nans
if flaginan
    H(INAN) = nan;
end

% prepare output
OUT = DEM;
OUT.Z = H;
OUT.name = ['local topography (' options.type ')'];

