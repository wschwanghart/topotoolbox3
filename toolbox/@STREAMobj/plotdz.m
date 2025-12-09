function h = plotdz(S,DEM,options)

%PLOTDZ Plot river profile (upstream distance versus elevation)
%
% Syntax
%
%     plotdz(S,DEM)
%     plotdz(S,DEM,pn,pv,...)
%     plotdz(S,z,pn,pv,...)
%     h = ...
%
% Description
%
%     Plot stream distance from the outlet versus elevation (or any other
%     variable).
%
% Input arguments
%
%     S      instance of STREAMobj
%     DEM    digital elevation model (GRIDobj)
%     z      node attribute list (as returned by various STREAMobj
%            methods, e.g. STREAMobj/streamorder, STREAMobj/gradient)
%
%     Parameter name/value pairs {default}
%
%     'annotation':  {[]} ix      
%     vector with linear indices of locations into the DEM. The cells
%     referenced by ix must be part of the stream network. Use snap2stream
%     to adjust locations. Annotation is achieved with vertical arrows.
%     Note that simple points to be plotted on the stream network is easier
%     using PPS/plotdz.
%
%     'annotationtext': cell array of strings
%     if annotated, a cell array of strings can be added to the vertical 
%     arrows
%
%     'distance': {S.distance}
%     node attribute list with custom distances (see STREAMobj/distance) or
%     STREAMobj (see function distance(S,S2)). 
%
%     'dunit': {'m'} 'km'
%     distance unit. plotdz assumes that distance is given in meter. 
%
%     'doffset': {0}
%     add an offset (scalar) to the distance from outlet to shift the 
%     x-axis of the plot.
%
%     'linewidth', 1
%     scalar that specifies the width of the line. 
%
%     'color': {'b'}
%     line colors as provided to plot. Alternatively, you can supply a node
%     attribute list (nal). The line will then have varying colors based on 
%     nal and h will be a surface object if 'colormethod' is 'surface'.
%
%     if 'color' is a node attribute list, then following parameter
%     name/values apply
%
%     'colormethod'    {'line'} or 'surface' or 'split'
%     lines with variable colors can be obtained by different methods.
%     Before introduction of the new graphics system in Matlab, using the
%     edges of surfaces was the standard hack. Since 2014b and the new
%     graphics system, lines have an undocumented edge property. Modifying
%     this edge property results in much smoother and more beautiful lines,
%     but may be bugged in newer versions. The option 'split' uses
%     splitbyattribute to split the line into multiple lines which are
%     then colored. In this case, the color is not related to the color
%     limits of the axis. The color limits of the axis will be adjusted to
%     the value range only if 'colorbar' is set to true.
%
%     'colormap'  {'parula'}
%     string that identifies a known colormap (e.g. 'jet','landcolor')
%
%     'colorbar'  {true} or false
%     true adds a colorbar. Note that colors created with the undocumented
%     method can not be changed afterwards by 'caxis'.
%
%     'cbarlabel' {''}
%     string to label colorbar
%     
% Output arguments
%
%     h     handle to the line handle. h will be a surface handle if color 
%           is set to a nal and colormethod is surface.
%
% Example 1
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'mex',true,'preprocess','carve');
%     S  = STREAMobj(FD,flowacc(FD)>1000);
%     S  = klargestconncomps(S);
%     plotdz(S,DEM)
%
% Example 2 (colored plot)
%     
%     A  = flowacc(FD);
%     c  = chitransform(S,A);
%     plotdz(S,DEM,'color',c,'cbarlabel','\chi [m]')
%
% Example 3 (adjust distance to be \chi)
%
%     plotdz(S,DEM,'distance',c)
%     xlabel('\chi [m]')  
%
% See also: STREAMobj, STREAMobj/plot, STREAMobj/smooth, PPS/plotdz
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 7. December, 2025

arguments
    S  STREAMobj
    DEM {mustBeGRIDobjOrNal(DEM,S)}
    options.parent = gca
    options.annotation = []
    options.color = getNextColor(gca)
    options.annotationtext = []
    options.distance  = S.distance
    options.dunit {mustBeMember(options.dunit,{'km' 'm'})} = 'm'
    options.doffset (1,1) = 0
    options.colormap = 'parula'
    options.linewidth (1,1) {mustBePositive} = 1
    options.colormethod {mustBeMember(options.colormethod,...
        {'line' 'surface' 'split'})} = 'line'
    options.colorbar (1,1) = true
    options.cbarlabel  = ''
    options.type {mustBeMember(options.type,...
        {'plot','area','stairs','stairsarea'})} = 'plot'
    options.EdgeColor = [.3 .3 .3]
    options.FaceColor = [.7 .7 .7]
    options.FaceAlpha (1,1) = 1
    options.EdgeAlpha (1,1) = 1
    options.ncolors (1,1) {mustBePositive,mustBeInteger} = 20
    options.BaseValue = []
    
end

ax   = options.parent;

% Get elevation nal
zz = ezgetnal(S,DEM);

% get dynamic properties of S
order    = S.orderednanlist;

% distance
if isempty(options.distance)
    dist = S.distance;
else
    if isa(options.distance,'STREAMobj')
        dist = distance(S,options.distance);
    elseif ischar(options.distance)
        dist = S.distance;
        dist = dist./max(dist);
    else
        dist = options.distance;
    end
end

% Convert horizontal distance to km?
switch lower(options.dunit)
    case 'km'
        dist = dist/1000;
end
        
% apply distance offset
dist = dist + options.doffset;

I     = ~isnan(order);
d     = nan(size(order));
d(I)  = dist(order(I));
z     = nan(size(order));
z(I)  = zz(order(I));

% Plot with uniform color
if ~isnal(S,options.color) && ~isa(options.color,'GRIDobj')
    switch options.type
        case 'plot'
            ht = plot(ax,d,z,'-','Color',options.color,...
                'LineWidth',options.linewidth);
        case 'stairs'
            ht = stairs(ax,d,z,'-','Color',options.color,...
                'LineWidth',options.linewidth);
        case 'area'
            if isempty(options.BaseValue)
                basevalue = min(z);
            else
                basevalue = options.BaseValue;
            end
            ht = area(ax,d,z,'EdgeColor',options.EdgeColor,...
                'FaceColor',options.FaceColor,...
                'FaceAlpha',options.FaceAlpha,...
                'EdgeAlpha',options.EdgeAlpha,...
                'BaseValue',basevalue);
        case 'stairsarea'
            if isempty(options.BaseValue)
                basevalue = min(z);
            else
                basevalue = options.BaseValue;
            end
            [xb,yb] = stairs(ax,d,z);
            ht = area(ax,flipud(xb),flipud(yb),'EdgeColor',options.EdgeColor,...
                'FaceColor',options.FaceColor,...
                'FaceAlpha',options.FaceAlpha,...
                'EdgeAlpha',options.EdgeAlpha,...
                'BaseValue',basevalue);
    end
    
else
    % Make sure that the colors come as node-attribute list  
    options.color = ezgetnal(S,options.color);
    switch lower(options.colormethod)
        case 'line'
            %% Plotting colored lines using undocumented Edges property
            % see here:
            % http://undocumentedmatlab.com/blog/plot-line-transparency-and-color-gradient
            
            ht    = plot(ax,d,z,'-');
            c     = zeros(size(order,1),3);
            minc  = min(+options.color);
            maxc  = max(+options.color);
            
            cmap    = colormap(ax,options.colormap)*255;
            cmapix  = linspace(minc,maxc,size(cmap,1));
            
            col   = interp1(cmapix,cmap,+options.color);
            c(I,:) = col(order(I),:);
            c     = c';
            c     = [c;zeros(1,size(c,2))+200];
            c     = uint8(c);
            % nans must be removed
            c     = c(:,I);
            
            ht.LineWidth = options.linewidth;
            % seems that the line must be first drawn ...
            drawnow
            % ... to be colored.
            set(ht.Edge, 'ColorBinding','interpolated', 'ColorData',c);
            
            if options.colorbar
                cc = colorbar(ax);
                clim([minc maxc]);
            end
            
            
        case 'surface'
            %% Plotting colored lines using surface
            
            dummy = z*0;
            c     = nan(size(order));
            c(I)  = +options.color(order(I));
            colormap(ax,options.colormap)
            ht = surface([d d],[z z],[dummy dummy],[c c],...
                'facecolor','none',...
                'edgecolor','flat',...
                'linewidth',options.linewidth,...
                'parent',ax);
            if options.colorbar
                cc = colorbar(ax);
            end
        otherwise
            [CS,zagg] = splitbyattribute(S,options.color,options.ncolors);
            [clr,lims] = num2rgb(zagg,options.colormap);
            if ishold(ax)
                keephold = true;
            else
                keephold = false;
            end
            hold(ax,'on');
            ht = cellfun(@(Ssub,col) ...
                plotdz(Ssub,DEM,'distance',nal2nal(Ssub,S,dist),...
                'color',col,'parent',ax,...
                'linewidth',options.linewidth),...
                CS(:),num2cell(clr,2));
            if options.colorbar
                colormap(ax,options.colormap);          
                hcb = colorbar(ax);
                if ~isempty(options.cbarlabel)
                    hcb.Label.String = options.cbarlabel;
                end
                clim(ax,lims)
            end

            if ~keephold
                hold(ax,'off')
            end
     
    end
    if options.colorbar && ~isempty(options.cbarlabel)
        cc.Label.String = options.cbarlabel;
    end
end

xlabel(['Distance upstream [' lower(options.dunit) ']'])
ylabel('Elevation [m]')

%% Annotation
if ~isempty(options.annotation)
    ix = options.annotation;
    hold on
    [Lia,Locb] = ismember(ix,S.IXgrid);


    if any(~Lia)
        error('TopoToolbox:STREAMobj',...
            'Some of the annotations are not located on the stream network')
    end
    
    annd = dist(Locb);
    
    if isa(DEM,'GRIDobj')
        annz = DEM.Z(S.IXgrid(Locb));
    else
        annz = zz(Locb);
    end
    
    if ~isempty(options.annotationtext)
        c = options.annotationtext;
        addtext = true;
    else
        addtext = false;
    end
    
    
    for r = 1:numel(ix)
        if addtext
            annotext = [c{r} '\newline\downarrow'];
        else
            annotext = '\downarrow ';
        end
        
        text(ax,'Position',[annd(r), annz(r)],...
             'String', annotext,...
             'VerticalAlignment','bottom',...
             'FontWeight','bold');
    end
    hold off
end


if nargout == 1
    h = ht;
end
end


function clr = getNextColor(ax)
        currentIndex = mod(ax.ColorOrderIndex - 1, size(ax.ColorOrder,1)) + 1;
        clr = ax.ColorOrder(currentIndex, :);
end