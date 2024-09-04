function h = plot(SW,options)
%PLOT Plot instance of SWATHobj
%
% Syntax
%
%     plot(S)
%     plot(S,...)
%     h = ...;
%
% Description
%
%     PLOT creates a map-view plot of a SWATHobj. The plot can be color-
%     coded by the z-values, usually elevation, for which either matlab's
%     'plot' or 'scatter' function can be used. This needs to be specified
%     in the parameters (see below) and affects the speed of the plotting,
%     which can be slow if the SWATHobj contains many data points.
%
% Input arguments
%
%     S     instance of SWATHobj
%
%     Parameter name/value pairs
%
%     'trace'   {true}, false
%     determines if the central trace of the SWATHobj is plotted
%
%     'outline'   {true}, false
%     determines if the outline of the SWATHobj is plotted
%
%     'points'   {true}, false
%     determines if the data points of the SWATHobj are plotted
%
%     'left'   {true}, false
%     determines if the left half (as seen from the direction of the
%     central line) of the SWATHobj is plotted
%
%     'right'   {true}, false
%     determines if the right half (as seen from the direction of the
%     central line) of the SWATHobj is plotted
%
%     'legend'   {true}, false
%     determines if a legend is plotted
%
%     'labeldist'   {[]}, scalar, 1 x m vector
%     adds labels to specific points along the profile. By default, no
%     labels are shown. If a scalar is provided, the profile is divided
%     into segments of equal lengths. If a vector is provided, only the
%     points at the specified distances are labeled.
%
%     'plotmode'   {'plot'}, 'scatter'
%     switches between Matlab's 'plot' and 'scatter' function for making
%     the plot. The scatter command is usually faster when 'colorz' is set
%     to true
%
%     'colorz'   true, {false}
%     optionally colors the data points by their z-values
%
%     'colorrange' 1x2 scalar vector
%     determines the minimum and maximum values of the colormap range.
%     Default values are the minimum and maxmum values of all data points
%     in the SWATHobj.
%
%     'colormap'   {'jet'}, one of Matlab's builtin colormaps
%     determines the colors used for plotting the data points if 'colorz'
%     is set to true. Note that if plotting on top of, e.g., a DEM produced
%     by imagesc(DEM), the colorbar of the DEM will change accordingly.
%
%     'colormode'   {'normal'}, 'inverse'
%     determines how the colors change with values. Default is 'normal',
%     i.e., low (high) z-values correspond to blue (red) colors in the
%     'jet' colormap.
%
%     'colorbar', {true}, false
%     optionally includes a colorbar; only if 'colorz' is true
%
%     'markersize'   {2}, scalar
%     size of the circular markers used to plot the location of the data
%     points
%
%
% Output arguments (optional)
%
%     h    handle object of the plotted features. If more than one feature
%          (trace, outline, points) is plotted, h contains more than one
%          handle
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     SW = SWATHobj(DEM,'dx',200,'dy',200);
%     G = gradient8(DEM,'degree');
%     SWG = mapswath(SW,G);
%     figure, plot(SWG,'colorz',true,'plotmode','scatter',...
%         'colorrange',[0 30],'colormap','hot','markersize',10,'legend',false);
%
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
% Date: May, 2015

arguments
    SW SWATHobj
    options.trace (1,1) = true
    options.tracecolor = 'k'
    options.outline (1,1) = true
    options.outlinecolor = 'k'
    options.points (1,1) = false
    options.left (1,1) = true
    options.right (1,1) = true
    options.legend (1,1) = true
    options.labeldist {mustBeNumeric} = []
    options.plotmode {mustBeMember(options.plotmode,{'plot','scatter'})} = 'plot'
    options.colorz (1,1) = false
    options.colormap = 'turbo'
    options.colorrange (1,2) = [-inf inf]
    options.colormode {mustBeMember(options.colormode,{'normal','inverse'})} = 'normal'
    options.colorbar (1,1) = true
    options.markersize (1,1) {mustBeNumeric,mustBePositive} = 2
end

labeldist  = options.labeldist;
colmap     = options.colormap;
colrange   = options.colorrange;
colmode    = options.colormode;
markersize = options.markersize;

% check
if ~options.left && ~options.right
    warning('Empty SWATHobj: nothing to plot')
    return;
end

% Create color codes
if (options.colorz)
    cmap = colormap(colmap)';
    if strcmp(colmode,'inverse'); cmap = flipud(cmap); end
    minz = colrange(1);
    maxz = colrange(2);
    if isinf(minz); minz = min(min(SW.Z)); end
    if isinf(maxz); maxz = max(max(SW.Z)); end
    CX = linspace(minz,maxz,length(cmap));
end

% Limit by side
if ~(options.left)
    ny = ceil(length(SW.disty)/2)+1;
    SW.X = SW.X(ny:end,:);
    SW.Y = SW.Y(ny:end,:);
    SW.Z = SW.Z(ny:end,:);
    
elseif ~(options.right)
    ny = floor(length(SW.disty)/2);
    SW.X = SW.X(1:ny,:);
    SW.Y = SW.Y(1:ny,:);
    SW.Z = SW.Z(1:ny,:);
end

% Plot SWATHobj data points
ct = 1;
if (options.points)    
    ix = find(~isnan(SW.Z));
    if~isempty(ix)
        xp = SW.X(ix);
        yp = SW.Y(ix);
        
        if (options.colorz)
            
            % Color-code z values
            zp = SW.Z(ix);
            
            switch options.plotmode
                case 'plot'
                    [N,edges,bins] = histcounts(zp,20);
                    bincenters = (edges(1:end-1)+edges(2:end))/2;

                    clr = num2rgb(bincenters,cmap');

                    for r = 1:numel(N)

                        I = bins == r;
                        hh = plot(xp(I),yp(I),'.',...
                            'Color',clr(r,:),'MarkerSize',markersize);
                        hold on
                        if r == 1
                            ht(ct) = hh;
                        end
                    end

                case 'scatter'
                    % Interpolate colors
                    mfcol = [interp1(CX,cmap(1,:),zp),...
                            interp1(CX,cmap(2,:),zp),interp1(CX,cmap(3,:),zp)];
                    ht(ct) = scatter(xp,yp,markersize,mfcol);
                    hold on
            end
            
        else
            
            % No color-coding
            mfcol = [0 0 0];
            switch options.plotmode
                case 'plot'
                    ht(ct) = plot(xp,yp,'.','color',mfcol,'Markersize',markersize); 
                    hold on
                case 'scatter'
                    ht(ct) = scatter(xp,yp,markersize,mfcol,'filled');
                    hold on
            end
        end
    end
    legstr(ct) = {'Points'};
    ct = ct+1;
end

% Trace of SWATHobj
if (options.trace)
    ht(ct) = plot(SW.xy0(:,1),SW.xy0(:,2),'o','Markersize',2,...
        'Color',options.tracecolor); 
    hold on
    legstr(ct) = {'Original trace'};
    ct = ct+1;
    ht(ct) = plot(SW.xy(:,1),SW.xy(:,2),'--','Markersize',2,...
        'Color',options.tracecolor);
    legstr(ct) = {'Resamp./interp. trace'};
    ct = ct+1;
end

% Outline of SWATHobj
if (options.outline)
    IM = ~isnan(SW.Z);
    B = bwboundaries(IM,4);
    for k = 1 : length(B)
        ix = sub2ind(size(IM),B{k}(:,1),B{k}(:,2));
        x_outline = SW.X(ix);
        y_outline = SW.Y(ix);
        ht(ct) = plot(x_outline,y_outline,'-','Color',options.outlinecolor);
    end
    legstr(ct) = {'Outline'};
end

% Label the distance
if ~isempty(labeldist)
    xp = SW.xy(:,1);
    yp = SW.xy(:,2);
    dp = SW.distx;
    if isscalar(labeldist)
        lp = (0:labeldist:max(dp))';
    elseif isvector(labeldist)
        lp = labeldist;
    end
    % interpolate label positions
    lx = interp1(dp,xp,lp);
    ly = interp1(dp,yp,lp);
    % add labels
    textstr = cellstr(num2str(lp));
    hold on, plot(lx,ly,'kv');%'wv')
    text(lx+2*SW.dx,ly,textstr,'color',[0 0 0]); % offset=2xSW.dx
end


drawnow
% Legend
if (options.legend)
    legend(ht,legstr);
end

% Colorbar
if (options.colorz) && (options.colorbar)
    hc = colorbar;
    yt = get(hc,'YTick');
    set(hc,'YTickLabel',(yt.*(maxz-minz)+minz)');
end

% Pass handles
if nargout>0
    h = ht;
end

end