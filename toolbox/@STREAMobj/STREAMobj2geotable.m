function GT = STREAMobj2geotable(S,varargin)

%STREAMobj2geotable Convert STREAMobj to a geotable
%
% Syntax
%
%     GT = STREAMobj2geotable(S)
%     GT = STREAMobj2shape(S,'seglength',seglength,'type',type)
%     GT = STREAMobj2shape(...,'attributes',{attname data aggfun ...})
%
% Description
%
%     STREAMobj2geotable takes a STREAMobj and converts it to a geotable. A
%     geotable (geospatial table) is a type of table designed to store and
%     manage geospatial data. It is specifically used for handling data
%     that includes geographic coordinates and associated attributes.
%     Geotables are part of MATLAB's Mapping Toolbox. They store coordinate
%     reference systems, and can be exported to shapefiles. Their use
%     should be preferred over mapstructs (see STREAMobj2mapstruct).
%
% Input parameters
%
%     S       STREAMobj
%     
%     Parameter name/value pairs
%
%     'seglength'   length of line features
%     'type'        'geo' or 'map'. Determines whether the function returns
%                   a geoshape or mapshape object.
%     'attributes'  cell array with attribute data (see STREAMobj2mapstruct
%                   for details). These data will be aggregated to feature
%                   properties. 
%
% Output arguments
%
%     GT      geotable
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     A = flowacc(FD);
%     c = chitransform(S,A);
% 
%     GT = STREAMobj2geotable(S,'type','geo',...
%             'attributes',{...
%             'z' DEM @mean ...
%             'chi' c @mean});
%
%     geoplot(GT,ColorVariable="chi");
% 
% See also: STREAMobj2mapstruct
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 13. June, 2024

d = S.distance;
defaultseglength = max(max(d)/20,5*S.cellsize);

% Check inputs
p = inputParser;
addParameter(p,'type','map',...
    @(x) ischar(validatestring(x,{'geo' 'map'},'STREAMobj2geotable')));
addParameter(p,'seglength',defaultseglength,@(x) isscalar(x) && x>S.cellsize);
addParameter(p,'attributes',{})
addParameter(p,'proj',projcrs())
parse(p,varargin{:});

% Calculate mapping structure
GT = STREAMobj2mapstruct(S,'seglength',p.Results.seglength,...
    'attributes',p.Results.attributes);

if isProjected(S)
    proj = S.georef.ProjectedCRS;
    t    = 'planar';
elseif isGeographic(S)
    proj = S.georef.GeographicCRS;
    t    = 'geographic';
else
    proj = [];
    t    = 'planar';
end

GT = mapstruct2geotable(GT,'CoordinateReferenceSystem',proj,...
    'coordinateSystemType',t);
