function STREAMobj2kml(S,filename,options)

%STREAMobj2kml Convert STREAMobj to kml (Google Earth)
%
% Syntax
%
%     STREAMobj2kml(S,filename)
%     STREAMobj2kml(S,filename,'pn','pv',...)
%     k = STREAMobj2kml(...)
%
% Description
%
%     STREAMobj2kml converts a stream network stored in a STREAMobj to kml
%     (or kmz) that can be opened in Google Earth and other GIS software.
%
%     STREAMobj2kml enables plotting of stream networks in Google Earth 
%     that are colored by an attribute (e.g. elevation, gradient, ...).
%
% Input arguments
%
%     S             STREAMobj
%     filename      filename including extension 
%
%     Parameter name/value pairs
%
%     'seglength'      Line segment length into which the stream network is
%                      split
%     'attribute'      GRIDobj or node-attribute list that will be used as
%                      attribute for each line segment. Default is
%                      S.distance. May not contain nans or infs.
%     'aggfun'         Anonymous function that will be used to aggregate
%                      node-attribute lists to segment attributes. By
%                      default, this is @mean.
%     'colormap'       string/char, anonymous function or colormap matrix
%                      used for plotting.
%     'nrcolors'       Number of colorvalues. Default is 20.
%     'linewidth'      Linewidth. Default is 2.
%     'name'           String/char. Name of the kmz file as it appears in
%                      Google Earth.
%
% Output arguments
%
%     k        kml object (see kml toolbox)
%
% Example: Plot a chimap in Google Earth
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     A = flowacc(FD);
%     c = chitransform(S,A);
%     [~,zb] = zerobaselevel(S,DEM);
%     c = c+zb;    
%     STREAMobj2kml(S,'test.kml','attribute',c)
%     % Open streamnet.kmz in Google Earth
%
% See also: STREAMobj2mapstruct, STREAMobj2shape, wmplot, kml
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 3. September, 2024

arguments
    S STREAMobj
    filename
    options.seglength (1,1) {mustBePositive} = inf
    options.attribute = []
    options.aggfun = @mean
    options.attributename = 'Attr'
    options.color  = ttclr('river')
    options.colormap = 'turbo'
    options.nrcolors (1,1) {mustBePositive,mustBeNumeric} = 20
    options.name = 'Stream network'
    options.linewidth (1,1) {mustBePositive,mustBeNumeric} = 2
end

if isempty(options.attribute)
    GT = STREAMobj2geotable(S,'type','geo','seglength',options.seglength);
    clr = options.color;
else
    if isinf(options.seglength)
        d = S.distance;
        seglength = max(max(d)/20,5*S.cellsize);
    else
        seglength = options.seglength;
    end

    GT = STREAMobj2geotable(S,'type','geo','seglength',seglength,...
        'attributes',...
        {options.attributename options.attribute options.aggfun});
    
    clr = num2rgb(GT.(options.attributename),options.colormap);
end

kmlwrite(filename,GT,"Color",clr,'LineWidth',options.linewidth,...
    'Name',options.name)
