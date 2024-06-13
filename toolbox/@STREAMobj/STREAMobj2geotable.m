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
GT = rmfield(GT,'Geometry');

% Make sure that vectors in GT are row vectors
fun = @(x) x(:)';
T = cellfun(fun,{GT.X},'UniformOutput',false);
[GT.X] = T{:};
T = cellfun(fun,{GT.Y},'UniformOutput',false);
[GT.Y] = T{:};


% Project to different coordinate system
if isProjected(S)
    switch lower(p.Results.type)
        case 'geo'
            for r = 1:numel(GT)
                [GT(r).Lat,GT(r).Lon] = projinv(S.georef.ProjectedCRS,...
                    GT(r).X,GT(r).Y);
            end
            GT = rmfield(GT,{'X','Y'});
            proj = geocrs(4326);
        otherwise
            proj = S.georef.ProjectedCRS;
    end
elseif isGeographic(S)
    switch lower(p.Results.type)
        case 'geo'

            for r = 1:numel(GT)
                GT(r).Lat = GT(r).Y;
                GT(r).Lon = GT(r).X;
            end
            GT = rmfield(GT,{'X','Y'});
            proj = geocrs(4326);

        otherwise

            for r = 1:numel(GT)
                [GT(r).X,GT(r).Y] = projfwd(p.Results.proj,...
                    GT(r).Y,GT(r).X);
            end            
            proj = p.Results.proj;

    end
else
    % Keep X and Y as they are
    proj = [];
end

% Convert to table
GT = struct2table(GT);

% Convert coordinates to *lineshapes
GT = table2geotable(GT,"GeometryType","line","CoordinateReferenceSystem",proj);

% Remove the coordinate variables
switch lower(p.Results.type)
    case 'geo'
    GT = removevars(GT,{'Lat','Lon'});
    case 'map'
    GT = removevars(GT,{'X','Y'});
end