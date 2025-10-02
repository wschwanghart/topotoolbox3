function [DEM,MASK] = widenstream(S,DEM,inp,method)

%WIDENSTREAM level elevations adjacent to the stream network
%
% Syntax
%
%     DEMw = widenstream(S,DEM,nrpx)
%     DEMw = widenstream(S,DEM,widthnal)
%     DEMw = widenstream(S,DEM,xyw,method)
%     DEMw = widenstream(S,DEM,'4conn')
%     DEMw = widenstream(S,DEM,'4conn',method)
%     [DEMw,MASK] = ...
%
% Description
%
%     This function sets elevations in a digital elevation model adjacent
%     to the stream network to the same elevations as the nearest stream
%     pixel. The function thus allows to level the valley bottom or stream
%     bed. Latter may be frequently needed when working with
%     high-resolution topographic data (e.g. LiDAR) and numerical,
%     hydrodynamic models (e.g. LISFLOOD FP).
%
%     widenstream(S,DEM,nrpx) levels elevations in DEM within a distance of 
%     nrpx pixels around the stream network S (STREAMobj).
%
%     widenstream(S,DEM,'4conn') ensures a downward descending flowpath
%     along cardinal stream pixels. This will ensure that flood models that
%     route along cardinal (and not diagonal) neighbors only will not
%     impound where the stream network S connects diagonal neighbors.
%     Diagonal neighbors can be connected by two cardinal neighbors. The
%     algorithm chooses and lowers the cardinal neighbor with the lower 
%     elevation while the other remains unchanged.
%
%     widenstream(S,DEM,'4conn',method) will adjust the elevation of
%     cardinal neighbors connecting diagonal neighbors using different
%     methods. 'linear' will take the average of the upstream and
%     downstream neighbors. 'previous' will adjust the elevation to that of
%     the upstream neighbor. 'next' will adjust the elevation to that of
%     the downstream neighbor. Default is 'linear'.
%
%     widenstream(S,DEM,xyw,method) levels elevations in DEM by the width
%     in mapunits measured at the locations xy. The locations and width are
%     must be provided as the nx3 matrix xyw where the first two columns
%     contain x and y coordinates, respectively, and the third column is
%     the width. The function snaps the measured locations to the 
%     stream network S and interpolates between them using
%     one-dimensional interpolation (interp1) with the method (all methods
%     accepted by the function interp1). 
%
%     *Note that widenstream(S,DEM,xyw,method) only works with a STREAMobj
%     of a single river reach, e.g. their must be only one channel head.*
%
% Input arguments
%
%     S        STREAMobj     
%     DEM      Digital Elevation Model (GRIDobj)
%     nrpx     number of pixels (half-widths)
%     widthnal node-attribute list containing the width for each stream
%              network node
%     xyw      n-by-3 matrix with x and y coordinates and river width.
%     method   interpolation method ('linear','spline',...). See interp1 
%              for further options. In case if option '4conn' is chosen,
%              method can be 'linear', 'previous', or 'next'.
%
% Output arguments
%
%     DEMw     Carved DEM.
%     MASK     Stream mask
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     S = trunk(S);
%     DEM = imposemin(S,DEM);
%     IX = [694708  553319  370738  262825  139925  1801]';
%     w  = [100 150 200 150 50 200]';
%     [x,y] = ind2coord(DEM,IX);
%     xyw = [x y w];
%     DEMw = widenstream(S,DEM,xyw,'linear');
%     imageschs(DEMw)
%
%
% See also: interp1, STREAMobj/imposemin 
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 2. October, 2025

arguments
    S  STREAMobj
    DEM GRIDobj
    inp 
    method {mustBeMember(method,{'linear','nearest','next','previous','spline','pchip','cubic'})} = 'linear'
end

% If users choose '4conn'
if (ischar(inp) || isstring(inp)) 
    if ~strcmpi(inp,"4conn")
        error('Unknown method for widening stream.')
    end

    switch method
        case 'linear'
            m = 1;
        case 'previous'
            m = 2;
        case 'next'
            m = 3;
        otherwise
            error('Method does not apply for 4conn.')
    end

    [row,col] = ind2sub(S.size,S.IXgrid);
    % Identify downstream cardinal neighbors
    I = isapprox(hypot(S.x(S.ix)-S.x(S.ixc),S.y(S.ix)-S.y(S.ixc)),...
                       S.cellsize);

    SG = STREAMobj2GRIDobj(S);

    for r = 1:numel(S.ix)
        
        if I(r)
            % Downstream neighbor is cardinal, proceed.
            continue
        end
    
        % Downstream neighbor is diagonal. There are two pixels that are
        % candidates to be lowered. 
        rix = row(S.ix(r));
        cix = col(S.ix(r));
        rixc = row(S.ixc(r));
        cixc = col(S.ixc(r));

        subs1   = [rix cixc];
        subs2   = [rixc cix];
        
        % Linear index of cardinal candidates
        ixd     = sub2ind(S.size,[subs1(1);subs2(1)],[subs1(2);subs2(2)]);

        % Which of the two has a lower elevation?
        [~,ixlowest] = min(DEM.Z(ixd));

        % What value should we assign to the pixel to be lowered
        % 1. Elevation of upstream pixel
        if m == 2
            DEM.Z(ixd(ixlowest)) = DEM.Z(S.IXgrid(S.ix(r)));
        
        % 2. Elevation of downstream pixel
        elseif m == 3
            DEM.Z(ixd(ixlowest)) = DEM.Z(S.IXgrid(S.ixc(r)));
        % 3. Average of upstream and downstream pixel
        else
        DEM.Z(ixd(ixlowest)) = (DEM.Z(S.IXgrid(S.ix(r))) + ...
                                 DEM.Z(S.IXgrid(S.ixc(r))))/2;
        end
        
        % Add pixel to the stream mask
        SG.Z(ixd(ixlowest)) = true;
    end

    if nargout == 2
        MASK = SG;
    end
    
    return
end

% if users choose the other options, we continue here
I = false(DEM.size);
I(S.IXgrid) = true;

[D,L] = bwdist(I,'e');

if isscalar(inp)
    I = D<=inp;
elseif isnal(S,inp)
    distnal = inp/DEM.cellsize;
    HG = zeros(DEM.size);
    HG(S.IXgrid) = distnal;
    I = HG(L) >= D;
else
    if numel(streampoi(S,'Channelheads','ix')) > 1
        error('TopoToolbox:widenstream',...
            ['widenstream works only with a river reach, i.e., there\n'...
             'must not be more than one channel head.'])
    end

    xyw = inp;
    xyw(:,3) = xyw(:,3)/DEM.cellsize/2;
    
    [~,~,IX] = snap2stream(S,xyw(:,1),xyw(:,2));
    d  = distance(S,'max_from_ch');
    [~,locb] = ismember(IX,S.IXgrid);
    
    nrObsPix = accumarray(locb,1);
    if any(nrObsPix>=2)
        warning('multiple observations per pixel, taking averages');
    end
    w    = accumarray(locb,xyw(:,3),size(S.IXgrid),@mean,nan);
    
    I    = isnan(w);
    w(I) = interp1(d(~I),w(~I),d(I),method,'extrap');
    W    = zeros(DEM.size);
    W(S.IXgrid) = w;
    
    I  = D<=W(L);
end
    
DEM.Z(I) = DEM.Z(L(I));

if nargout == 2
    MASK = GRIDobj(DEM,'logical');
    MASK.Z = I;
end
