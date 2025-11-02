function [GT,removed] = clipshape(GT,x,y)

%CLIPSHAPE Clip geotable to xy-limits or polygon
%
% Syntax
%
%     GTc = clipshape(GT,shape)
%     GTc = clipshape(GT,xlimits,ylimits)
%     GTc = clipshape(GT,latlims,lonlims)
%     [GTc,removed] = clipshape(GT,...)
%
% Description
%
%     clipshape clips a geotable with lines, polygons or points to the
%     extent in the clipper polygon shape (a map- or geopolyshape object).
%     The clipper polygon can also be defined by a rectangle defined by
%     xlimits and ylimits or latitude and longitude limits. 
%
% Input arguments
%
%     GT        Geotable (it is assumed that the geotable consists of only 
%               one type of geometries, e.g. polygons, lines or points).
%     shape     clipper polygon (must have same coordinate reference system 
%               as GT)
%     xlimits   limits of x-coordinate (two element vector)
%     ylimits   limits of y-coordinate (two element vector)
%     latlims   limits of latitudes (two element vector)
%     lonlims   limits of longitudes (two element vector)
%
% Output arguments
%
%     GTc       clipped geotable
%     removed   logical vector ([size(GT,1) 1]) with true elements
%               indicating features that were outside the clipper polygon 
%               that have been removed.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     D   = drainagebasins(FD);
%     GT  = GRIDobj2geotable(D);
%     ext = getextent(pad(DEM,-200));    
%     [GTc,removed] = clipshape(GT,ext(1:2),ext(3:4));
%     mapshow(GT)
%     hold on
%     mapshow(GTc,'FaceColor','b')
%
%     S = STREAMobj(FD,'minarea',1000);
%     P = PPS(S,'rpois',0.001);
%     [GTs,GTp] = as(P,'geotable');
%     GTsc = clipshape(GTs,union(GTc.Shape));
%     GTpc = clipshape(GTp,union(GTc.Shape));
%     mapshow(GTsc,'color','w')
%     mapshow(GTpc,'color','w')
% 
% See also: mapclip, geoclip
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 2. November, 2025

arguments
    GT
    x
    y = []
end

if isstruct(GT)
    GT = mapstruct2geotable(GT);
end

[~,isproj] = parseCRS(GT);

if isempty(y)

    if isproj
        shapes = mapclip(GT.Shape,x);
    else
        shapes = geoclip(GT.Shape,x);
    end
else
    if isproj
        shapes = mapclip(GT.Shape,x,y);
    else
        shapes = geoclip(GT.Shape,x,y);
    end
end

HasVertices = false(size(GT,1),1);
for r = 1:size(GT,1)
    switch GT.Shape(1).Geometry
        case 'point'
            HasVertices(r) = shapes(r).NumPoints > 0;
        case 'line'
            HasVertices(r) = shapes(r).NumParts > 0;
        otherwise
            HasVertices(r) = shapes(r).NumRegions > 0;
    end
end

GT.Shape = shapes;
GT = GT(HasVertices,:);

if nargout == 2
    removed = ~HasVertices;
end

