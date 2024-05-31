function varargout = getoutline(DEM,options)

%GETOUTLINE Get outline of GRIDobj
%
% Syntax
%
%     MS = getoutline(DEM)
%     [x,y] = getoutline(DEM)
%     ... = getoutline(DEM,'pn',pv')
%
% Description
%
%     getoutline returns the outline of a GRIDobj. By default, getoutline
%     returns the coordinate vectors of the DEM edges. By setting shownans
%     = true, you can get the outline around the valid (non-nan) data in
%     the DEM.
%
% Input arguments
%
%     DEM         GRIDobj
%     
%     Parameter name/value pairs
%     
%     'shownans'    {false} or true
%     'output'      {mappolyshape}, mapshape, mapstruct (only applicable if
%                   called with one output argument.
%
% Output arguments
%
%     MS     mapping structure that can be displayed using mapshow or
%            exported with shapewrite (alternatively, see option 'output'.
%     x,y    coordinate vectors that can be used to plot the extent
%            rectangle (plot(x,y))
%
% Example 1
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEM.Z(DEM.Z>1500) = nan;
%     MS  = getoutline(DEM,'shownans',true,'output','mappolyshape');
%     % Plot in geographic axes or mapaxes object
%     geoplot(MS)
%
% Example 2
%
%     DEM = readexample('taalvolcano');
%     I = identifyflats(DEM);
%     I.Z = bwareaopen(I.Z,20);
%     DEM.Z(I.Z) = nan;
%     
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 30. May, 2024


arguments (Input)
    DEM  GRIDobj
    options.shownans = 0
    options.output = "mappolyshape" 
    options.simplify = true
end

% Check if there are nans
I = isnan(DEM);
% It's only required to account for nans if there are nans in the grid
nnan = options.shownans && any(I);

if ~nnan
    csh   = DEM.cellsize/2;
    [x,y] = getcoordinates(DEM);
    maxx = max(x)+csh/2;
    minx = min(x)-csh/2;
    
    maxy = max(y)+csh/2;
    miny = min(y)-csh/2;
    
    x = [minx minx maxx maxx minx];
    y = [miny maxy maxy miny miny];

    xy = [x(:) y(:)];

else
    I = ~I;
    [B,~,~,A] = bwboundaries(I.Z,8,"holes","TraceStyle","pixeledge",...
        "CoordinateOrder","xy");
    G   = digraph(A');
    deg = indegree(G);
    d   = distances(G,find(deg == 0));
    isenclosed = mod(d-1,2) == 0;

    % flip boundaries of enclosed boundaries
    for r = 1:numel(B)
        if isenclosed(r)
            B{r} = flipud(B{r});
        end
    end

    % get coordinates
    wf = DEM.wf;
    B = cellfun(@(cr)[(wf*[cr-1 ones(size(cr,1),1)]')';[nan nan]],B,...
        'UniformOutput',false);
    
    xy = vertcat(B{:});
end

% simplify lines
if options.simplify
    xy = dpsimplify(xy,100*eps);
end
x = xy(:,1);
y = xy(:,2);

% Without output arguments, the outline will be plotted
if nargout == 0
    if ~wm
        plot(x,y,'LineWidth',3,'Color',[.8 .8 .8]);
        hh = ishold;
        hold on
        plot(x,y,'k--','LineWidth',1);
        if ~hh
            hold off;
        end
    else
        try
        [lat,lon] = projinv(DEM.georef.ProjectedCRS,x,y);
        catch
            lat = y;
            lon = x;
        end
        wmline(lat,lon)

    end
        
elseif nargout == 1
    % One output argument, create mapping structure
    % split at nans
    switch options.output
        case 'geotable'
            MS = mappolyshape(x,y);
            MS.ProjectedCRS = DEM.georef.ProjectedCRS;
            MS = table(MS,1,'VariableNames',{'Shape','ID'});
        case 'mappolyshape'
            MS = mappolyshape(x,y);
            MS.ProjectedCRS = DEM.georef.ProjectedCRS;
        case 'mapshape'
            MS = mapshape(x,y,'Geometry','polygon');
        case 'mapstruct'

            if ~isnan(x(end))
                x = [x;nan];
                y = [y;nan];
            end
        
            I = isnan(x);
            nrlines = nnz(I);
            ix = find(I)';
            ixs = [1; ix(1:end-1)'+1];
            ixe = ix-1;
            for r = 1:nrlines
        
                MS(r).Geometry = 'Polygon';
                MS(r).X = x(ixs(r):ixe(r));
                MS(r).Y = y(ixs(r):ixe(r));
                MS(r).ID = r;
        
            end

        case 'geopolyshape'
            [lat,lon] = projinv(DEM.georef.ProjectedCRS,x,y);
            MS = geopolyshape(lat,lon);
            MS.GeographicCRS = geocrs(4326);

        case 'geoshape'
            [lat,lon] = projinv(DEM.georef.ProjectedCRS,x,y);
            MS = geoshape(lat,lon,'Geometry','polygon');
        case ''
    end
    varargout{1} = MS;

elseif nargout == 2
    % Two outputs, return x and y coordinates
    
    varargout{1} = x;
    varargout{2} = y;
end

end
    
       