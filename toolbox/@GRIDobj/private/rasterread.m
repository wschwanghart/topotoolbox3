function [dem,wf,header] = rasterread(file)

% read/import ESRI ascii grid
%
% Syntax
%
%     [dem,wf] = rasterread
%     [dem,wf] = rasterread('file')
%     [dem,wf,header] = rasterread('file')
%
% Description
%
%     read ASCII raster exported from ESRI ArcGIS. 'file' is a string 
%     indicating location and name of the file, e.g. 'srtm.txt'. When
%     called with no input arguments or a folder path, rasterread opens a
%     dialog box for retrieving the file.
% 
%     the header must contain following rows (in arbitrary order)
%     _________________________
%     ncols         201
%     nrows         451
%     xllcorner     329995.5
%     yllcorner     5297997.5
%     cellsize      5
%     NODATA_value  -9999
%     _________________________
% 
%     dem is a matrix with ncols columns and nrows rows. X and Y are 
%     coordinate matrices and provide the spatial reference for the dem.
%     header is a 6-x-1 structur array that contains the information
%     of the header of the file.
%
%
% See also: RASTERWRITE, DLMREAD, DLMWRITE
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 25 November, 2010



% read header data
persistent lastfolder


if nargin == 0
    
    filedialogbox = true;
    file = lastfolder;
    
elseif nargin == 1
    
    % check if only a folder is supplied
    if isfolder(file) || isempty(file)
        filedialogbox = true;
    else
        filedialogbox = false;
        if exist(file,'file') ~= 2
            error('file or folder does not exist')
            
        end
    end
    
end
    
if filedialogbox
    FilterSpec  = {'*.txt';'*.asc'};
    DialogTitle = 'Select ESRI ASCII grid';
    [FileName,PathName] = uigetfile(FilterSpec,DialogTitle,file);
    
    if FileName == 0
        error('TopoToolbox:incorrectinput',...
              'no file was selected')
    end
    
    file = [PathName FileName];
    lastfolder = PathName;
end

if ischar(file)
else
    error('file is a string indicating the ASCII file name and position')
end



fid=fopen(file,'r');
% loop through header

header = struct('ncols',[],...
                'nrows',[],...
				'xllcorner',[],...
				'yllcorner',[],...
				'cellsize',[],...
				'nodata',[]);
				
names   = fieldnames(header);
nrnames = numel(names);

try
fseek(fid,0,'bof');
for r = 1:nrnames 
    headertext = fgetl(fid);
    [headertext, headernum] = strtok(headertext,' ');
    I = cellfun(@(x,y) strcmpi(x(1:4),y(1:4)),names,repmat({headertext},nrnames,1));
	header.(names{I}) = str2double(headernum);
end
catch ME1
    error('header can not be read')
end


% read raster data
dem = fscanf(fid,'%lg',[header.ncols header.nrows]);
fclose(fid);
dem(dem==header.nodata) = NaN;

dem = dem';


% create world-file
wf = zeros(2,3);
wf(1,1) = header.cellsize;
wf(2,2) = -header.cellsize;
wf(1,3) = header.xllcorner + header.cellsize/2;
wf(2,3) = header.yllcorner + header.nrows*header.cellsize - header.cellsize/2;


