function out = ttupdate(options)

%TTUPDATE Automatically update TopoToolbox addon
%
% Syntax
%
%     result = ttupdate(verbose = true, nolibtt = false)
%
% Description
%
%     This function attempts to automatically update TopoToolbox with the
%     latest release available on github. If the function fails, use 
%     clear all
%     before running this function. 
%
%     Note that the function will remove your current version. Any
%     changes to toolbox files will be lost.
%
% Input arguments
%
%     Parameter name/value pairs
%
%     'verbose'   {true} or false. If true, information on the progress 
%                 will be displayed in the command window
%     'nolibtt'   true or {false}. If true, a MATLAB-only version will be
%                 installed without libtopotoolbox.
% 
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 23. November, 2025

arguments
    options.verbose (1,1) = true
    options.nolibtt (1,1) = false
end

verbose = options.verbose;

% Github information
owner = "topotoolbox";
repo  = "topotoolbox3";
url   = sprintf("https://api.github.com/repos/%s/%s/releases/latest", ...
                owner, repo);

if verbose
    disp(string(datetime("now")) + " - Get download URL")
end
% Read API to get the latest 
opts = weboptions("ContentType","json","Timeout",10);
json  = webread(url,opts);

if verbose
    disp(string(datetime("now")) + " - Get operating system")
end

if ~options.nolibtt
    archstr = computer('arch');
end

if verbose && ~options.nolibtt
    disp(string(datetime("now")) + " - OS is " + string(archstr))
end

assets  = json.assets;
if ~options.nolibtt
    I = contains({assets.name},archstr);
else
    I = contains({assets.name},"nolibtt");
end

if verbose
    disp(string(datetime("now")) + ...
        " - Download (" + string(round(assets(I).size/1e6,1)) + " MB).");
end
directory   = tempdir;
toolboxfile = fullfile(directory,assets(I).name);
file = websave(toolboxfile,...
    assets(I).browser_download_url);

if verbose
    disp(string(datetime("now")) + " - Download successful to " + string(file))
end

if verbose
    disp(string(datetime("now")) + " - Deinstall old version.")
end

% Find addons
addons = matlab.addons.toolbox.installedToolboxes;

% Uninstall previous version of toolbox
for k = 1:numel(addons)
    if addons(k).Name == "TopoToolbox"
        matlab.addons.toolbox.uninstallToolbox(addons(k));
    end
end

% Install new version
if verbose
    disp(string(datetime("now")) + " - Install new version.")
end
out = matlab.addons.toolbox.installToolbox(toolboxfile);

if verbose
    disp(string(datetime("now")) + " - Delete temporary toolbox file.")
end
delete(toolboxfile)

if verbose
    disp(string(datetime("now")) + " - Installation successful.")
end

