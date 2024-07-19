function ht = wmplot(S,options)

%WMPLOT plot stream network in the webmap browser
%
% Syntax
%
%     h = wmplot(S)
%     h = wmplot(S,'color',c)
%
% Description
%
%     wmplot plots the stream network S in MATLAB's webmap browser. This
%     requires the Mapping Toolbox and S must have a valid georeferencing.
%
%     Note that wmplot displays colored stream networks, but that the
%     values are classified rather than continuous.
%
% Input arguments
%
%     S    STREAMobj
%
% Output arguments
%
%     h    handle to wmline object
%     
%     Parameter name/value pairs
%
%     'color'      three element rgb vector or node attribute list
%     'colormap'   char (default is 'jet')
%     'nrcolors'   number of colors (only applicable, if color is a node
%                  attribute list)
%
% Example 1
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     wmplot(S)
%
% Example 2: Plot ksn in a webmap.
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     A = flowacc(FD);
%     k = ksn(S,DEM,A,0.45,20);
%     wmplot(S,'color',k)
%
% Remark
%
%     Another way to quickly plot a stream object in a webmap is:
%
%     [lat,lon] = STREAMobj2latlon(S);
%     h = wmline(lat,lon,'OverlayName','Stream network','color','k');
%
% See also: STREAMobj, STREAMobj/plot, STREAMobj/plotc,
%           STREAMobj/splitbyattribute
%           
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 18. July, 2024

arguments
    S  STREAMobj
    options.color = 'k'
    options.colormap = 'turbo'
    options.nrcolors  = 20
end

if ~isnal(S,options.color) && ~isa(options.color,'GRIDobj')
    % STREAMobj to lat lon
    if isGeographic(S)
        [lon,lat] = STREAMobj2XY(S);
    elseif isProjected(S)
        [lat,lon] = STREAMobj2latlon(S);
    else
        error("S has an undefined CRS.")
    end
    
    h = wmline(lat,lon,'OverlayName','Stream network','Color',options.color);
    
    if nargout == 1
        ht = h;
    end
    
else

    [CS,zagg] = splitbyattribute(S,options.color,options.nrcolors);

    clr = num2rgb(zagg,options.colormap);


    counter = 0;
    h = cell(1,numel(CS));
    for r = 1:numel(CS)
        counter = counter + 1;
        h{counter} = wmplot(CS{r},'color',clr(r,:));
    end
    if nargout == 1
        ht = horzcat(h{:});
    end
end
