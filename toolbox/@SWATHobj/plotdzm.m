function h = plotdzm(SW,M,options)
%PLOTDZM create color-coded distance-elevation plot from SWATHobj and GRIDobj
%
% Syntax
%
%     plotdzm(SW,M)
%     plotdzm(SW,M,'pn','pv',...)
%     h = ...;
%
% Description
%
%     PLOTDZM creates a profile-view plot of a SWATHobj, showing the
%     z-values as a function of distance along the profile of the SWATHobj
%     and color-coded by the z-values in the GRIDobj M.
%
%     For display the function can use either Matlab's 'plot' or 'scatter'
%     function or create an image. This needs to be specified in the
%     parameters (see below) and affects the speed of the plotting, which
%     can be slow if the SWATHobj contains many data points. In this case,
%     creating an image is fastest.
%
% Input arguments
%
%     S     instance of SWATHobj
%     M     instance of GRIDobj
%
%     Parameter name/value pairs
%
%     'left'   {true}, false
%     determines if the left half (as seen when looking along the direction
%     of the central line) of the SWATHobj is plotted.
%
%     'right'   {true}, false
%     determines if the right half (as seen when looking along the
%     direction of the central line) of the SWATHobj is plotted.
%
%     'direction'   {'x'}, 'y'
%     plot data either along the x-direction (default) of the SWATHobj or
%     along the y-direction (across the SWATHobj).
%
%     'distadjust'   {0}, scalar
%     allows shifting the x-axis by a scalar value. This can be useful when
%     alligning SWATHobj's that were obtained along a certain reach of
%     a drainage network, with the STREAMobj that it was derived from.
%
%     'distance'    {}, 1 x 2 vector
%     choose minimum and maximum distance along the SWATHobj for plotting.
%
%     'zmode'   {'absolute'}, 'relative'
%     assigns elevations along the SWATHobj either in absolute values
%     (default), or relative to the minimum value across each line in the
%     SWATHobj (relative). Note that this can be inaccurate if the points
%     across the SWATHobj miss to sample the lowest elevation present.
%
%     'sortm'   {'descend'}, 'ascend', 'none'
%     determines if and how the values within each component of a SWATHobj
%     are sorted before being plotted. This parameter affects the
%     visibility of specific populations of the data points if the markers
%     overlap in the plot.
%
%     'colormap'   {'gray'}, one of Matlab's builtin colormaps
%     determines the colors used for plotting the data points if 'colorz'
%     is set to true. Note that if plotting on top of, e.g., a DEM produced
%     by imagesc(DEM), the colorbar of the DEM will change accoringly.
%
%     'colormode'   {'normal'}, 'inverse'
%     determines how the colors change with values. Default is 'normal',
%     i.e., low (high) m-values would correspond to blue (red) colors in
%     the 'jet' colormap.
%
%     'colorrange'   1x2 vector
%     determines the minimum and maximum values of the colormap range.
%     Default values are the minimum and maxmum values of the entire
%     SWATHobj.
%
%     'colorbar'   true, {false}
%     adds a colorbar to the plot. If plotmode is set to 'plot' or 'scatter'
%     the default colorbar that is created with the corresponding button of
%     a matlab figure does not apply to the plotted data. If plotmode is
%     set to 'image' this feature works.
%
%     'plotmode'   'plot', 'scatter', {'image'}
%     defines how the data points are plotted. If the number of data is
%     rather small 'plot' and 'scatter' may produce a figure that can be
%     exported as editable vector graphics. This feature breaks down at
%     larger numbers of data points and the time required for building the
%     figure can get very(!) long. In such cases it more useful to
%     interpolate an 'image'.
%
%     'markersize'   {2}, scalar
%     size of the circular markers used to plot the location of the data
%     points. Only applicable if 'plotmode' set to 'plot' or 'scatter'.
%
%     'xstretch' and 'ystretch'    {1 and 1}, scalar
%     if 'plotmode' is set to 'image', the x- and y-values can be stretched
%     by the corresponding factors to increase the resolution of the image.
%   
%     'backgrd'     'high','low',{'none'}
%     if set to 'high' or 'low', the background is given the highest or
%     lowest value in the data range, respectively. This helps to control
%     the appearance of the background when the plotmode is 'image'.
%
% 
% Output arguments (optional)
%
%     h    handle object of the plotted features. If more than one feature
%          (trace, outline, points) is plotted, h contains more than one
%          handle
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     FA = flowacc(FD);
%     S = STREAMobj(FD,FA>1e6/FA.cellsize/FA.cellsize);
%     S = klargestconncomps(S);
%     S = trunk(S);
%     [x,y] = STREAMobj2XY(S);
%     ix = ~isnan(x);
%     SW = SWATHobj(DEM,x(ix),y(ix),'dx',100,'dy',100,'width',4000,'smooth',2000);
%     G = gradient8(DEM,'degree');
%     figure,plotdzm(SW,G,'colorrange',[0 35],'ystretch',0.1), colorbar
%     title('Hillslope angles within 2 km of the river')
%
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
%         Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 15. September, 2024

arguments
    SW    SWATHobj
    M     GRIDobj
    options.left (1,1) = true
    options.right (1,1) = true;
    options.direction {mustBeTextScalar,mustBeMember(options.direction,{'x','y'})} = 'x';
    options.distadjust (1,1) {mustBeNumeric} = 0
    options.distance {mustBeNumeric} = []
    options.zmode {mustBeTextScalar,mustBeMember(options.zmode,{'absolute','relative'})} = 'absolute'
    options.colormap = 'pink'
    options.colormode {mustBeTextScalar,mustBeMember(options.colormode,{'normal','inverse'})} = 'normal'
    options.colorrange (1,2) = [-inf inf]
    options.colorbar (1,1) = false
    options.colorbarlabel  = ''
    options.sortm {mustBeTextScalar,mustBeMember(options.sortm,{'descend','ascend','none'})} = 'descend'
    options.plotmode  {mustBeTextScalar,mustBeMember(options.plotmode,{'plot','scatter','image'})} = 'image'
    options.xstretch (1,1) {mustBeNumeric,mustBePositive} = 1
    options.ystretch (1,1) {mustBeNumeric,mustBePositive} = 1
    options.markersize (1,1) {mustBeNumeric,mustBePositive} = 4
    options.backgrd  {mustBeTextScalar,mustBeMember(options.backgrd,{'none','low','high'})} = 'none'
end

% parameters
left       = options.left;
right      = options.right;
direction  = options.direction;
distadjust = options.distadjust;
distance   = options.distance;
zmode      = options.zmode;
colmap     = options.colormap;
colmode    = options.colormode;
colrange   = options.colorrange;
colbar     = options.colorbar;
cbarlabel  = options.colorbarlabel;
sortm      = options.sortm;
plotmode   = options.plotmode;
markersize = options.markersize;
xstretch   = options.xstretch;
ystretch   = options.ystretch;
backgrd    = options.backgrd;


warning('off')
[SWG] = mapswath(SW,M);
warning('on')

if ~left && ~right
    warning('Empty SWATHobj: nothing to plot')
    return;
end

minm = colrange(1);
maxm = colrange(2);
if isinf(minm); minm = min(min(SWG.Z)); end
if isinf(maxm); maxm = max(max(SWG.Z)); end

if ~strcmp(plotmode,'image')
    hf = figure('visible','off');
    cmap = colormap(colmap)';
    % cb = colorbar;
    % cdata = get(get(cb,'Children'),'Cdata');
    close(hf)
    if strcmp(colmode,'inverse'); cmap = fliplr(cmap); end
    CX = linspace(minm,maxm,length(cmap));
end


if ~left
    ny = ceil(length(SW.disty)/2)+1;
    SW.Z = SW.Z(ny:end,:);
    SWG.Z = SWG.Z(ny:end,:);
elseif ~right
    ny = floor(length(SW.disty)/2);
    SW.Z = SW.Z(1:ny,:);
    SWG.Z = SWG.Z(1:ny,:);
end

if ~isempty(distance)
    ix = SW.distx<distance(1) | SW.distx>distance(2);
    SW.Z(:,ix) = nan;
    SWG.Z(:,ix) = nan;
end

if strcmp(zmode,'relative')
    minv = min(SW.Z,[],1);
    SW.Z = SW.Z-repmat(minv,length(SW.disty),1);
end

[zr,zc] = size(SW.Z);
switch direction
    case 'x'
        d = repmat(SW.distx',zr,1);
        swdx = SW.dx;
    case 'y'
        d = repmat(SW.disty,1,zc);
        swdx = SW.dy;
end
z = reshape(SW.Z,numel(SW.Z),1);
s = reshape(SWG.Z,numel(SWG.Z),1);
d = reshape(d,numel(d),1);

d = d(~isnan(z));
s = s(~isnan(z));
z = z(~isnan(z));

% Sort
if ismember(sortm,{'ascend','descend'})
    [s,ix] = sort(s,sortm);
    z = z(ix);
    d = d(ix);
end

% Interpolate colors
if ~strcmp(plotmode,'image')
    mfcol = [interp1(CX,cmap(1,:),s),...
        interp1(CX,cmap(2,:),s),interp1(CX,cmap(3,:),s)];
end

switch plotmode
    case 'plot'
        for k = 1 : length(z)
            plot(d(k)+distadjust,z(k),'s','color',mfcol(k,:),'Markerfacecolor',mfcol(k,:),'Markersize',markersize), hold on
        end
        
    case 'scatter'
        hout = scatter(d+distadjust,z,markersize,mfcol,'square','filled'); hold on
        
    case 'image'
        x = d+distadjust;
        y = z;
        v = s;
        
        % apply colorrange thresholds
        x = x(v>=minm);
        y = y(v>=minm);
        v = v(v>=minm);
        v = min(v,maxm);
        y = round(y.*ystretch); % stretch y-values
        x = round(x.*xstretch); % stretch x-values
        
        F = TriScatteredInterp(x,y,v);
        
        xy = unique([x,y],'rows');
        qv = F(xy(:,1),xy(:,2));
        
        qx = floor(min(x)):swdx:ceil(max(x));
        %             qx = floor(min(x)):ceil(max(x));
        qy = floor(min(y)):ceil(max(y));
        IM = GRIDobj(qx./swdx,qy,nan(length(qy),length(qx)));
        %             IM = GRIDobj(qx,qy,zeros(length(qy),length(qx)));
        ix = coord2ind(IM,xy(:,1)./swdx,xy(:,2));
        %             ix = coord2ind(IM,xy(:,1),xy(:,2));
        IM.Z(ix(~isnan(ix))) = qv(~isnan(ix));
        %             IM.Z(ix) = qv;
        
        %             SE = strel('square',3)./9;
        %             IM.Z = imdilate(IM.Z,SE);
        %             SE = ones(5)./25;
        %             IM.Z = conv2(IM.Z,SE);
        
        switch backgrd
            case 'high'
                IM.Z(isnan(IM.Z)) = maxm;
            case 'low'
                IM.Z(isnan(IM.Z)) = minm;
        end
        
        
        hout = imagesc(qx./xstretch,flipud(qy'./ystretch),IM.Z,[minm maxm]);
        set(gca,'YDir','normal');
        hc = colormap(colmap);
        if strcmp(colmode,'inverse')
            colormap(flipud(hc));
        end
        axis normal
end



drawnow
xlabel(sprintf('Distance along profile (%s)',SW.xyunit))
ylabel(sprintf('Z (%s)',SW.zunit))

if colbar && ~strcmp(plotmode,'image')
    colormap(colmap);
    hc = colorbar;
    if strcmp(colmode,'inverse')
        set(get(hc,'Children'),'CData',(64:-1:1)');
    end
    yt = get(hc,'YTick');
    set(hc,'YTickLabel',(yt.*(maxm-minm)+minm)');
    if not(isempty(cbarlabel))
        hc.Label.String = cbarlabel;
        hc.Label.FontSize = 10;
        hc.FontSize = 10;
    end
elseif colbar && strcmp(plotmode,'image')
    hc = colorbar;
end

if nargout==1
    h = unique(hout);
    if exist('hc','var')
        h(end+1) = hc;
    end
end

