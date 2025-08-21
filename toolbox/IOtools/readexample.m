function DEM = readexample(example,options)

%READEXAMPLE read DEM from TopoToolbox DEM github repository
%
% Syntax
%
%     DEM = readexample(example)
%
% Description
%
%     readexample reads DEMs from the TopoToolbox DEM repository.
%     Currently available examples are
%
%     'greenriver'
%     'kedarnath'
%     'kunashiri'
%     'perfectworld'
%     'taalvolcano'
%     'taiwan'
%     'tibet'
%
%     The function requires internet connection.
%
% Input arguments
%
%     example    string of example name (e.g. 'tibet')
%
% Output arguments
%
%     DEM        example data. This is usually a GRIDobj. However, in some
%                cases, DEM is a structure array with variable data.
%     
% See also: GRIDobj, websave, readopentopo
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 1. August, 2025

arguments
    example 
    options.filename = [tempname '.tif']
    options.deletefile (1,1) = true
    options.verbose (1,1) = true
end

example = lower(example);

% create file to which the data will be saved
f = fullfile(options.filename);

switch example
    case 'taiwan'
        url = 'https://github.com/TopoToolbox/DEMs/raw/master/taiwan.tif';
        istif = true;
    case 'tibet'
        url = 'https://github.com/TopoToolbox/DEMs/raw/master/tibet.tif';
        istif = true;
    case 'taalvolcano'
        url = 'https://github.com/TopoToolbox/DEMs/raw/master/taalvolcano.tif';
        istif = true;
    case 'kunashiri'
        url = 'https://github.com/TopoToolbox/DEMs/raw/master/kunashiri.tif';
        istif = true;
    case 'perfectworld'
        url = 'https://github.com/TopoToolbox/DEMs/raw/master/perfectworld.tif';
        istif = true;
    case 'kedarnath'
        url = 'https://github.com/TopoToolbox/DEMs/raw/master/kedarnath.tif';
        istif = true;
    case 'bigtujunga'
        url = 'https://github.com/TopoToolbox/DEMs/raw/master/bigtujunga.tif';
        istif = true;
    case 'greenriver'
        url = 'https://github.com/TopoToolbox/DEMs/raw/master/greenriver.tif';
        istif = true;
    otherwise 
        error('There is no such example file.')
end
      
% save to drive
webopts = weboptions('Timeout',100000);

% Download with websave
if options.verbose
    disp([char(datetime("now")) ' -- Downloading...'])
end
outfile = websave(f,url,webopts);
if options.verbose
    disp([char(datetime("now")) ' -- Download finished...'])
end

% Read grid or mat file
if istif
    DEM = GRIDobj(f);
    DEM.name = example;
else
    DEM = load(f); 
end

if options.deletefile
    delete(f);
    if options.verbose
        disp([char(datetime("now")) ' -- Temporary file deleted'])
    end
end