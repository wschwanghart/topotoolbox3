function filename = openRasterDialog(folder)

%Opens a dialog window to select a file

persistent lastfolder

if nargin == 0
    if exist('lastfolder','var')
        folder = lastfolder;
    else        
        folder = pwd;
    end
end

% for a txt or tiff file as input
FilterSpec  = {'*.txt;*.asc;*.tif;*.tiff','supported file types (*.txt,*.asc,*.tif,*.tiff)';...
    '*.txt',   'ESRI ASCII grid (*.txt)';...
    '*.asc',   'ESRI ASCII grid (*.asc)';...
    '*.tif',   'GeoTiff (*.tif)';...
    '*.tiff',  'GeoTiff (*.tiff)';...
    '*.*',     'all files (*.*)'};

DialogTitle = 'Select ESRI ASCII grid or GeoTiff';
[FileName,PathName] = uigetfile(FilterSpec,DialogTitle,folder);

if FileName == 0
    error('TopoToolbox:incorrectinput',...
        'no file was selected')
end

filename = fullfile(PathName, FileName);
lastfolder = PathName;
