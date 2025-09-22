function GT = buffer(S,d,options)

%BUFFER Buffer polygon around stream network
%
% Syntax
%
%     GT = buffer(S,d)
%     GT = buffer(S,d,pn,pv,...)
%
% Description
%
%     Buffer creates a polygon around the stream network at the distance
%     d. By default, the function returns a geotable. 
%
% Input arguments
%
%     S      STREAMobj
%     d      buffer distance (must have the same units as the STREAMobj,
%            usually meters). d can also be a node-attribute list that
%            contains a buffer distance for each node in the network. To
%            run with reasonable times, the function classifies these
%            distances into 100 bins. Depending on the variability of
%            values in d, buffer distances might thus be slightly different
%            from values given in d. 
% 
%     Parameter name/value pairs
%
%     'output'      'geotable' or 'polyshape'. Default is 'geotable' and 
%                   requires the Mapping Toolbox.
%     'JointType'   Default is 'round'. See polybuffer for additional 
%                   options.
%
% Output arguments
%
%     GT     geotable or polyshape object
%
% Example 1: Uniform buffer distance
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S   = STREAMobj(FD,minarea = 1000);
%     S   = klargestconncomps(S,1);
%     GT  = buffer(S,100);
%     geoplot(GT)
%
% Example 2: Buffer distance as a function of upstream area
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S   = STREAMobj(FD,minarea = 1000);
%     S   = klargestconncomps(S,1);
%     A   = flowacc(FD);
%     d   = sqrt(getnal(S,A);
%     GT  = buffer(S,d);
%     geoplot(GT)
%
% See also: mappolyshape, polybuffer
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 12. September, 2025


arguments
    S  STREAMobj
    d {mustBePositive}
    options.output {mustBeMember(options.output,["geotable","polyshape"])} = "geotable"
    options.JointType = "round"
end

if isGeographic(S)
    error("The stream network must have a projected coordinate system.")
end

if ~isnal(S,d)
    [x,y] = STREAMobj2XY(S);
    polyout = polybuffer([x,y],"lines",d,"JointType",options.JointType);
else
    [CS,dagg] = splitbyattribute(S,d,100);
    [Cx,Cy] = cellfun(@(S) STREAMobj2XY(S),CS,'UniformOutput',false);
    Cpoly = cell(size(Cx));
    parfor r = 1:numel(Cx)
        Cpoly{r} = polybuffer([Cx{r},Cy{r}],"lines",dagg(r),...
            "JointType",options.JointType);
    end
    polyout = vertcat(Cpoly{:});
    
    polyout = union(polyout);
end


switch options.output
    case "geotable"

        nreg = polyout.NumRegions;
        R    = regions(polyout);
        
        % Preallocate geotable
        GT   = table(mappolyshape(cell(nreg,1),cell(nreg,1)),(1:nreg)','VariableNames',{'Shape' 'ID'});

        for k = 1:nreg
            % Convert vertices to mappolyshape
            mpg  = mappolyshape(R(k).Vertices(:,1),R(k).Vertices(:,2));  
            % Write mappolyshape in geotable
            GT.Shape(k) = mpg;
        end
        
        % Assign coordinate reference system
        GT.Shape.ProjectedCRS = parseCRS(S);

    otherwise
        GT = polyout;
end
end

