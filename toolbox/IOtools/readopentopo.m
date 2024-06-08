function DEM = readopentopo(options)

%READOPENTOPO Read global DEMs using the opentopography.org API
%
% Syntax
%
%     DEM = readopentopo(pn,pv,...)
%
% Description
%
%     readopentopo reads DEMs from opentopography.org using the API
%     described on:
%     http://www.opentopography.org/developers
%     Most of the DEMs come in geographic coordinates (WGS84) and should be
%     projected to a projected coordinate system (use reproject2utm or 
%     project) before analysis in TopoToolbox.     
%
%     NOTE: Starting on January 1st, 2022, an API authorization key will be 
%     required for this API. Users can request an API key via myOpenTopo in 
%     the OpenTopography portal (https://opentopography.org/developers).
%     See also the description below to learn how to make your API key
%     permanently available to readopentopo.
%
% Input arguments
%
%     Parameter name values
%     'interactive'    {true} or false. If true, readopentopo will open a
%                      GUI that enables interactive selection. If true,
%                      then any given extent options will be ignored.
%     'filename'       provide filename. By default, the function will save
%                      the DEM to a temporary file in the system's temporary 
%                      folder. The option 'deletefile' controls whether the
%                      file is kept on the hard drive.
%     'extent'         GRIDobj or four element vector with geographical 
%                      coordinates in the order [west east south north].
%                      If a GRIDobj is supplied, readopentopo uses the
%                      function GRIDobj/getextent to obtain the bounding
%                      box in geographical coordinates. If extent is set,
%                      then the following parameter names 'north',
%                      'south', ... are ignored.
%     'addmargin'      Expand the extent derived from 'extent',GRIDobj by a
%                      scalar value in °. Default is 0.01. The option is
%                      only applicable if extent is provided by a GRIDobj.
%     'north'          northern boundary in geographic coordinates (WGS84).
%                      The option is ignored if the option 'extent' is
%                      provided or if 'interactive', true.
%     'south'          southern boundary
%     'west'           western boundary
%     'east'           eastern boundary
%     'demtype'        The global raster dataset *
%                      {'SRTMGL3'}:       SRTM GL3 (90m) (default)
%                      'SRTMGL1':         SRTM GL1 (30m)  
%                      'SRTMGL1_E':       SRTM GL1 (Ellipsoidal)  
%                      'AW3D30':          ALOS World 3D 30m  
%                      'AW3D30_E':        ALOS World 3D (Ellipsoidal)
%                      'SRTM15Plus':      Global Bathymetry SRTM15+ V2.1  
%                                         (only mediterranean area so far)
%                      'NASADEM':         NASADEM Global DEM 
%                      'COP30':           Copernicus Global DSM 30m 
%                      'COP90':           Copernicus Global DSM 90m 
%                      'EU_DTM:           Continental Europe Digital 
%                                         Terrain Model 
%                      'GEDI_L3':         Global Ecosystem Dynamics 
%                                         Investigation 1x1 km DTM
%                      'GEBCOIceTopo':    Global Bathymetry (500 m)
%                      'GEBCOSubIceTopo': Global Bathymetry (500 m)
%                        
%                      * requires API Key (see option 'apikey').
%
%     'apikey'         char or string. Users can request an API key via 
%                      myOpenTopo in the OpenTopography portal. You can
%                      also create a text file in the folder IOtools named
%                      opentopography.apikey which must contain the API
%                      Key. If there is a file that contains the key, there 
%                      is no need to provide it here.
%     'verbose'        {true} or false. If true, then some information on
%                      the process is shown in the command window.
%     'checkrequestlimit' {true} or false. Opentopography implements
%                      request limits. If the chose extent exceeds the 
%                      request limit, readopentopo will issue an error. If
%                      this option is set to false, readopentopo will try
%                      to download.
%     'deletefile'     {true} or false. True, if file should be deleted
%                      after it was downloaded and added to the workspace.
% 
% Output arguments
%
%     DEM            Digital elevation model in geographic coordinates
%                    (GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEM2 = readopentopo('extent',DEM);
%     DEM2 = reproject2utm(DEM2,90);
%     imagesc(DEM2)
%     hold on
%     getoutline(DEM)
%     hold off
%
% See also: GRIDobj, websave, roipicker
%
% Reference: http://www.opentopography.org/developers
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 6. June, 2024

arguments
    options.demtype {mustBeTextScalar} = 'SRTMGL3'
    options.filename = [char(tempname) '.tif']
    options.interactive = false
    options.extent = []
    options.addmargin = 0.01
    options.north = 37.091337
    options.south = 36.738884
    options.west  = -120.168457
    options.east  = -119.465576
    options.deletefile = true
    options.verbose = true
    options.apikey  = ''
    options.checkrequestlimit = true
end

validdems = {'SRTMGL3','SRTMGL1','SRTMGL1_E',...
     'AW3D30','AW3D30_E','SRTM15Plus',...
     'NASADEM','COP30','COP90',...
     'EU_DTM','GEDI_L3','GEBCOSubIceTopo','GEBCOIceTopo'};

demtype = validatestring(options.demtype,...
    validdems,'readopentopo');

% Access global topographic datasets including SRTM GL3 (Global 90m), 
% GL1 (Global 30m), ALOS World 3D and SRTM15+ V2.1 (Global Bathymetry 500m). 
% Note: Requests are limited to 
% 125,000,000 km2 for SRTM15+ V2.1,'GEBCOSubIceTopo','GEBCOIceTopo', and 
% 4,050,000 km2 for SRTM GL3, COP90 and 
% 450,000 km2 for all other data.
requestlimits = [4.05e6, 0.45e6, 0.45e6, ...
                 0.45e6, 0.45e6, 125e6,...
                 0.45e6, 0.45e6, 4.05e6, ...
                 0.45e6, 50e7, 125e6, ...
                 125e6]; % km^2
requestlimit  = requestlimits(strcmp(demtype,validdems));

% API URL
url = 'https://portal.opentopography.org/API/globaldem?';

% create output file
f = fullfile(options.filename);

% check api

if isempty(options.apikey)
    % check whether file opentopography.apikey is available
    if exist('opentopography.apikey','file')
        fid = fopen('opentopography.apikey');
        apikey = textscan(fid,'%c');
        apikey = apikey{1}';
        % Remove trailing blanks, if there are any
        apikey = deblank(apikey);
    else
        error('Readopentopo requires an API Key. Please read the help.')
    end
else
    apikey = options.apikey;
end


% save to drive
woptions = weboptions('Timeout',100000);

% get extent
if ~isempty(options.extent)
    if isa(options.extent,'GRIDobj')
        ext = getextent(options.extent,true);
        west  = ext(1) - options.addmargin;
        east  = ext(2) + options.addmargin;
        south = ext(3) - options.addmargin;
        north = ext(4) + options.addmargin;
        
    elseif numel(options.extent) == 4
        west = options.extent(1);
        east = options.extent(2);
        south = options.extent(3);
        north = options.extent(4);
    else
        error('Unknown format of extent')
    end
else
    west = options.west;
    east = options.east;
    south = options.south;
    north = options.north;
end

% now we have an extent. Or did the user request interactively choosing
% the extent.
if any([isempty(west) isempty(east) isempty(south) isempty(north)]) || options.interactive
    if options.interactive 
        ext = roipicker('requestlimit',requestlimit);
        if isempty(ext)
            DEM = [];
            return
        end
        west = ext(2);
        east = ext(4);
        south = ext(1);
        north = ext(3);
    end
end

% Check request limit
a = areaint([south south north north],...
            [west east east west],...
            almanac('earth','radius','kilometers'));

if options.checkrequestlimit && (a > requestlimit)
    error('TopoToolbox:readopentopo',...
            ['Request limit (' num2str(requestlimit) ' km^2) exceeded.\n'...
             'Your extent is ' num2str(a,1) ' km^2. Choose a smaller area.' ])
end

if options.verbose

    disp('-------------------------------------')
    disp('readopentopo process:')
    disp(['DEM type: ' demtype])
    disp(['API url: ' url])
    disp(['Local file name: ' f])
    disp(['Area: ' num2str(a,2) ' sqkm'])
    disp('-------------------------------------')
    disp(['Starting download: ' char(datetime("now"))])
end

% Download with websave
if isempty(apikey)
    outfile = websave(f,url,'west',west,...
        'east',east,...
        'north',north,...
        'south',south,...
        'outputFormat', 'GTiff', ...
        'demtype', demtype, ...
        woptions);
else
    outfile = websave(f,url,'west',west,...
        'east',east,...
        'north',north,...
        'south',south,...
        'outputFormat', 'GTiff', ...
        'demtype', demtype, ...
        'API_Key',apikey,...
        woptions);
end

if options.verbose
    disp(['Download finished: ' char(datetime("now"))])
    disp(['Reading DEM: ' char(datetime("now"))])
end

try

    DEM      = GRIDobj(f);
    
    if ~isProjected(DEM)
        disp(' ')
        disp('The downloaded DEM is not in a projected coordinate system.')
        disp('Make sure to project the DEM using GRIDobj/project or')
        disp('GRIDobj/reproject2utm.')
        disp('  ')
    end
    
    DEM.name = demtype;
    if options.verbose
        disp(['DEM read: ' char(datetime("now"))])
    end
    
catch
    % Something went wrong. See whether we can derive some information.
    fid = fopen(outfile);
    in = textscan(fid,'%c');
    disp('Could not retrieve DEM. This is the message returned by opentopography API:')
    disp([in{1}]')
    disp('readopentopo returns empty array')
    fclose(fid);
    DEM = [];
end
    

if options.deletefile
    delete(f);
    if options.verbose
        disp('Temporary file deleted')
    end
end

if options.verbose
    disp(['Done: ' char(datetime("now"))])
    disp('-------------------------------------')
end
end
