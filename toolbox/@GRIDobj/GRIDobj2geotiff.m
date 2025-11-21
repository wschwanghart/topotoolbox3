function GRIDobj2geotiff(A,file)

%GRIDobj2geotiff Export an instance of GRIDobj to a geotiff file
%
% Syntax
%    
%     GRIDobj2geotiff(DEM)
%     GRIDobj2geotiff(DEM,filename)
%
% Description
%
%     GeoTIFF is a common image file format that stores coordinates and
%     projection information to be read by most GIS software.
%     GRIDobj2geotiff writes an instance of GRIDobj to a GeoTIFF file. 
%
%     GRIDobj2geotiff requires the function geotiffwrite available with 
%     the Mapping Toolbox. If geotiffwrite does not exist on the search
%     path, the function will write a standard tif together with a
%     '.tfw'-file (worldfile, http://en.wikipedia.org/wiki/World_file ) to
%     the disk. 
%
%     GRIDobj2geotiff(DEM) opens a dialogue box to save the GeoTIFF
%
%     GRIDobj2geotiff(DEM,filename) saves the DEM to the specified 
%     filename
%
% Input arguments
%
%     DEM        instance of GRIDobj
%     filename   absolute or relative path and filename
%
% See also: GRIDobj
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 14. November, 2025

arguments
    A  GRIDobj
    file = ''
end

% if only 1 argument, open file dialog box
if isempty(file)
    [FileName,PathName] = uiputfile({'*.tif'});
    if FileName == 0
        disp('     no output written to disk')
        return
    end
    file = [PathName FileName];
end

% try to use geotiffwrite, which comes with the Mapping Toolbox
try
    if isempty(A.georef)
        % There might be no georef entry. In this case, we will try to
        % write a geographic grid as no information about the projection is
        % available
        R = georasterref(A.wf,A.size);
        geotiffwrite(file,A.Z,R,"CoordRefSysCode","EPSG:4326");
    else
        % A.georef is not empty. If a CRS is available, the following code
        % will use geotiffwrite
        if isProjected(A) || isGeographic(A)

            GeoKeyDirectoryTag = createGeoKeyDirectoryTag(A);
            geotiffwrite(file,A.Z,A.georef,...
                "GeoKeyDirectoryTag",GeoKeyDirectoryTag ...
                );
        else
            geotiffwrite(file,A.Z,A.georef,"CoordRefSysCode","EPSG:4326");
        end

    end
catch ME
    warning('TopoToolbox:GRIDobj',...
        ['GRIDobj2geotiff is unable to write a geotiff. Either you don''t \n'...
         'have the mapping toolbox, or there was another issue with geotiffwrite. \n'...
         'GRIDobj2geotiff instead writes a tif-image together with a world \n'...
         'file (*.tfw) which contains data on spatial referencing of the \n' ...
         'image, yet which lacks information on the type of projection used.']); 
              
    % if geotiffwrite is not available or any other error occurs
    % a tif file will be written to the disk together with a worldfile
    % .tfw-file.
    [pathstr, name, ~] = fileparts(file);
    dlmwrite(fullfile(pathstr,[name '.tfw']),A.wf(:),'precision', '%.10f');
    A = A.Z;
     
    siz = size(A);
    cla = class(A);
    
    switch cla
        case 'double'
            BpS = 64;
            TSF = Tiff.SampleFormat.IEEEFP;
        case 'single'
            BpS = 32;
            TSF = Tiff.SampleFormat.IEEEFP;
        otherwise
            if islogical(A)
                A = uint32(A);
                cla = 'uint32';
            end
            BpS = round(log2(double(intmax(cla))));
            TSF = Tiff.SampleFormat.UInt;
            
    end
    
    t = Tiff(file,'w');
    tagstruct.ImageLength = siz(1);
    tagstruct.ImageWidth = siz(2);
    tagstruct.BitsPerSample = BpS;
    tagstruct.SampleFormat = TSF;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.RowsPerStrip = 16;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    tagstruct.Photometric = 0;
    
    t.setTag(tagstruct);
    t.write(A);
    t.close;
end

end


function GeoKeyDirectoryTag = createGeoKeyDirectoryTag(DEM)
% Attempts to create a GeoKeyDirectoryTag from the referencing information
% in the maprefcells object in R

isproj = isProjected(DEM);

% Projected or geographic?
if isproj
    GeoKeyDirectoryTag.GTModelTypeGeoKey = 1; % Projected coordinate system
    GeoKeyDirectoryTag.GTCitationGeoKey = DEM.georef.ProjectedCRS.Name;
    GeoKeyDirectoryTag.GeogCitationGeoKey = DEM.georef.ProjectedCRS.GeographicCRS.Name;
    wkt = wktstring(DEM.georef.ProjectedCRS);
    wkt = char(wkt);

    epsg = extractEPSG(wkt);

    % ix = strfind(wkt,'"EPSG",');
    % c  = wkt(ix(end):end);
    % c  = extractAfter(c,'"EPSG",');
    % c  = extractBefore(c,']');

    GeoKeyDirectoryTag.ProjectedCSTypeGeoKey = epsg(end); %str2double([hemisphere sprintf('%02d',str2double(zone(regexp(zone,'[0-9]'))))]);
    GeoKeyDirectoryTag.ProjLinearUnitsGeoKey = 9001; % Linear_Meter
else
    GeoKeyDirectoryTag.GTModelTypeGeoKey = 2;
    GeoKeyDirectoryTag.GTCitationGeoKey = DEM.georef.GeographicCRS.Name;
    GeoKeyDirectoryTag.GeogCitationGeoKey = DEM.georef.GeographicCRS.Name;
    % wkt = wktstring(DEM.georef.GeographicCRS);
end

% Reference = cells or postings?
switch lower(DEM.georef.RasterInterpretation)
    case 'cells'
        GeoKeyDirectoryTag.GTRasterTypeGeoKey = 1; % RasterPixelIsArea
    case 'postings'
        GeoKeyDirectoryTag.GTRasterTypeGeoKey = 0; % RasterPixelIsPoint
end

GeoKeyDirectoryTag.GeogAngularUnitsGeoKey = 9102; %Angular_Degree

end

function epsgCode = extractEPSG(wktString)
    % This function extracts the EPSG code from a WKT string.
    %
    % Input:
    % wktString - A string containing the WKT representation of the CRS.
    %
    % Output:
    % epsgCode - The EPSG code as a numeric value. Returns NaN if EPSG code is not found.

    % Regular expression to find the EPSG code in the WKT string
    pattern = 'ID\["EPSG",(\d+)\]';
    matches = regexp(wktString, pattern, 'tokens');
    if ~isempty(matches)
        epsgCode = cellfun(@str2double,matches);
    else
        epsgCode = nan;
    end
end 