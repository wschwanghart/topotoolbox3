function data = readopenalti(options)

%READOPENALTI Read ICESAT altimetry data using the openaltimetry.org API
%
% Syntax
%
%     data = readopenalti(pn,pv,...)
%
% Description
%
%     readopenalti reads altimetry data from openaltimetry.org using the 
%     API described on: https://openaltimetry.earthdatacloud.nasa.gov/data
%     The data comes in geographic coordinates (WGS84). EGM96 geoid heights
%     are added using the Mapping Toolbox function egm96geoid if available.
%
% Input arguments
%
%     Parameter name values
%     'filename'       provide filename. By default, the function will save
%                      the data to a temporary file in the system's temporary 
%                      folder. The option 'deletefile' controls whether the
%                      file is kept on the hard drive.
%     'date'           datetime scalar or vector indicating the days
%                      required. Default is the last 600 days.
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
%     'north'          northern boundary in geographic coordinates (WGS84)
%     'south'          southern boundary
%     'west'           western boundary
%     'east'           eastern boundary
%     'product'        'atl03' L2A Global Geolocated Photon Data
%                      'atl06' L3A Land Ice Height
%                      'atl07' L3A Sea Ice Height
%                      'atl08' L3A Land and Vegetation Height (default)
%                      'atl10' L3A Sea Ice Freeboard
%                      'atl12' L3A Ocean Surface Height
%                      'atl13' L3A Inland Water Surface Height, Version 1                     
%     'level3a'        {true} or false
%     'verbose'        {true} or false. If true, then some information on
%                      the process is shown in the command window
%     'deletefile'     {true} or false. True, if file should be deleted
%                      after it was downloaded and added to the workspace.
% 
% Output arguments
%
%     data     table with altimetry data
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     d = readopenalti('extent',DEM,'addmargin',0,...
%           'date',datetime(2020,05,1):datetime('today'));
%     % coordinates are returned as longitude and latitude
%     [d.x,d.y] = mfwdtran(DEM.georef.mstruct,d.latitude,d.longitude);
%     surf(DEM,'block',true);
%     colormap(landcolor)
%     camlight
%     axis off
%     hold on
%     % plot points with some offset to elevation so that points are
%     % plotted above the surface
%     plot3(d.x,d.y,d.h_te_best_fit+20,'.r')
%
% See also: GRIDobj, websave, readopentopo, egm96geoid
%
% Reference: https://openaltimetry.earthdatacloud.nasa.gov/data/openapi/swagger-ui/index.html
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 12. September, 2024

arguments
    options.filename = [tempname '.csv']
    options.interactive = false
    options.addmargin = 0.01 %°
    options.extent   = []
    options.north (1,1) {mustBeNumeric,mustBeInRange(options.north,-90,90)} = 37.091337
    options.south (1,1) {mustBeNumeric,mustBeInRange(options.south,-90,90)} = 36.738884
    options.west  (1,1) {mustBeNumeric,mustBeInRange(options.west,-180,180)} = -120.16845
    options.east  (1,1) {mustBeNumeric,mustBeInRange(options.east,-180,180)} = -119.46557
    options.date = datetime('today')-600 : datetime('today')
    options.product {mustBeTextScalar} = 'atl08'
    options.level3a (1,1)    = false
    options.deletefile (1,1) = true
    options.verbose (1,1)    = true
end

validproducts = {'atl03','atl10','atl12',...
                 'atl13','atl06','atl07',...
                 'atl08'};

product = validatestring(options.product,validproducts,'readopenalti');
urltracks = 'https://openaltimetry.earthdatacloud.nasa.gov/data/api/icesat2/getTracks?';
if ~options.level3a
    url = ['https://openaltimetry.earthdatacloud.nasa.gov/data/api/icesat2/' product '?'];
else
    url = 'https://openaltimetry.earthdatacloud.nasa.gov/data/api/icesat2/level3a?';
end

% create output file
f = fullfile(options.filename);
    
% save to drive
webopts = weboptions('Timeout',100000);

% get extent
if ~isempty(options.extent)
    if isa(options.extent,'GRIDobj')
        % getextent with second input arg true returns lat/lon extent
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
    
    ext = roipicker();
   
    west = ext(2);
    east = ext(4);
    south = ext(1);
    north = ext(3);
end
    
if options.verbose
    a = areaint([south south north north],...
                [west east east west],almanac('earth','radius','kilometers'));
    disp('-------------------------------------')
    disp('readopenalti process:')
    disp(['Product product: ' product])
    disp(['API url tracks: ' urltracks])
    disp(['API url product: ' url])
    disp(['Local file name: ' f])
    disp(['Area: ' num2str(a,2) ' sqkm'])
    disp('-------------------------------------')
    disp(['Starting download: ' char(datetime('now'))])
end

% get date
if isdatetime(options.date)
    dt = options.date;
    D  = string(dt,'yyyy-MM-dd');
else
    D = options.date;
    dt = datetime(D,'InputFormat','yyyy-MM-dd');
end

if options.verbose
    disp(['Get tracks: ' char(datetime('now'))])
end

% Download tracks with websave
tracktable = table([],[],'VariableNames',{'track','date'});
if options.verbose
    disp(['Identifying tracks: ' char(datetime('now'))])
end

for r=1:numel(dt)
	   textwaitbar(r,numel(dt),'Wait');
	   
	   [~] = websave(f,urltracks,...
              'date',D(r),...
              'minx',west,...
              'maxx',east,...
              'maxy',north,...
              'miny',south,...
              'client','TopoToolbox',...
              'outputFormat', 'csv', ...
              webopts);
    temptracks = readtable(f);
    temptracks.date = repmat(dt(r),size(temptracks,1),1);
    tracktable = [tracktable; temptracks]; %#ok<AGROW>
    delete(f);
	   
end


for r = 1:numel(dt)
    
end

totaltracksfound = size(tracktable,1);

if options.verbose
    disp([num2str(totaltracksfound) ' tracks found.'])
    disp(['Data download starts: ' char(datetime('now'))])
end

% Download data
if options.level3a
    
    counter = 1;
    trackIds = unique(tracktable.track);
    startDate = D(1);
    endDate = D(end);
    for r = 1:numel(trackIds)
        trackId = trackIds(r);
        if options.verbose
            disp(['Download trackId ' num2str(trackId) ': ' char(datetime('now'))])
        end
        try
            [~] = websave(f,url,...
                'product',product,...
                'startDate',startDate,...
                'endDate',endDate,...
                'minx',west,...
                'maxx',east,...
                'maxy',north,...
                'miny',south,...
                'trackId',trackId,...
                'client','TopoToolbox',...
                'outputFormat', 'csv', ...
                webopts);
            if counter == 1 || ~exist('data','var')
                data = readtable(f);
            else
                data = [data; readtable(f)]; %#ok<AGROW>
            end
            counter = counter + 1;
        catch
            disp(['Download failed: ' char(datetime('now'))])
        end
        delete(f);
    end
    
else
    
    counter = 1;
    for r = 1:totaltracksfound
        dd = string(tracktable.date(r),'yyyy-MM-dd');
        trackId = num2str(tracktable.track(r));
        if options.verbose
            disp(['Download ' char(dd) ',' num2str(trackId) ': ' char(datetime('now'))])
        end
        try
            [~] = websave(f,url,...
                'date',dd,...
                'minx',west,...
                'maxx',east,...
                'maxy',north,...
                'miny',south,...
                'trackId',trackId,...
                'client','TopoToolbox',...
                'outputFormat', 'csv', ...
                options);
            if counter == 1  || ~exist('data','var')
                data = readtable(f);
            else
                data = [data; readtable(f)]; %#ok<AGROW>
            end
            counter = counter + 1;
        catch
            disp(['Download failed: ' char(datetime('now'))])
        end
        delete(f);
    end
end

if options.verbose
    disp(['Total of ' num2str(size(data,1)) ' points downloaded: ' char(datetime('now'))])

end

if options.verbose
    disp(['Adding EGM96 heights: ' char(datetime('now'))]);
end

data.egm96geoid = egm96geoid(data.latitude,data.longitude);

if options.verbose
    disp(['Done: ' char(datetime('now'))])
    disp('-------------------------------------')
end
end


function textwaitbar(i, n, msg)
% A command line version of waitbar.
% Usage:
%   textwaitbar(i, n, msg)
% Input:
%   i   :   i-th iteration.
%   n   :   total iterations.
%   msg :   text message to print.
%
% Date      : 05/23/2019
% Author    : Xiaoxuan He   <hexxx937@umn.edu>
% Institute : University of Minnesota
%
% Previous percentage number.
persistent i_prev_prct;
% Current percentage number.
i_prct = floor(i ./ n * 100);
% Print message when counting starts.
if isempty(i_prev_prct) || i_prct < i_prev_prct
    i_prev_prct = 0;
    S_prev = getPrctStr(i_prev_prct);
    
    fprintf('%s: %s',msg, S_prev);
end
% Print updated percentage.
if i_prct ~= i_prev_prct
    S_prev = getPrctStr(i_prev_prct);
    fprintf(getBackspaceStr(numel(S_prev)));
    
    S = getPrctStr(i_prct);
    fprintf('%s', S);
    
    i_prev_prct = i_prct;
end
% Clear percentage variable.
if i_prct == 100
    fprintf(' Done.\n');
    clear i_prev_prct;
end
end
function S = getPrctStr(prct)
S = sprintf('%d%%  %s',prct,getDotStr(prct));
if prct < 10
    S = ['  ',S];
elseif prct < 100
    S = [' ',S];
end
end
function S = getDotStr(prct)
S = repmat(' ',1,10);
S(1:floor(prct/10)) = '.';
S = ['[',S,']'];
end
function S = getBackspaceStr(N)
S = repmat('\b',1,N);
end

