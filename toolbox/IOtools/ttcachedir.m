function cachedir = ttcachedir(options)
%TTCACHEDIR Return the platform-specific cache directory for TopoToolbox 3
%
% Description
%
%     The default storage location for downloaded data is ttcachedir. Use
%     ttclearcache to delete files from the cache. The default storage
%     location is
%
%     Windows:
%         %LOCALAPPDATA%\topotoolbox
%
%     macOS:
%         $HOME/Library/Caches/topotoolbox
%
%     Linux:
%         $XDG_CACHE_HOME/topotoolbox
%
%     If this directory does not exist, ttcachedir will try to create it.
%     If another location should be used for the cache, supply the cachedir
%     argument.
%
% Input arguments
%
%     Parameter name/value pairs
%
%     'cachedir'   char or string. If supplied, ttcachedir will use this
%     directory instead of the default one.
%
% Author: William Kearney (william.kearney[at]uni-potsdam.de)
% Date: 30 March 2026

arguments
    options.cachedir = ''
end

if isempty(options.cachedir)

    if ispc()
        if ~isenv("LOCALAPPDATA")
            error("Please set environment variable LOCALAPPDATA or supply " + ...
                "cachedir argument to ttcachedir");
        end
        cachedir = getenv("LOCALAPPDATA");
    elseif isunix()
        if ismac()
            if isenv("HOME")
                cachedir = fullfile(getenv("HOME"),"Library","Caches");
            else
                error("Please set environment variable HOME or supply " + ...
                    "cachedir argument to ttcachedir");
            end
        else
            cachedir = getenv("XDG_CACHE_HOME");
            if isempty(cachedir)
                if isenv("HOME")
                    cachedir = fullfile(getenv("HOME"), ".cache");
                else
                    error("Please set environment variable HOME or supply " + ...
                        "cachedir argument to ttcachedir");
                end
            end
        end
    end

    cachedir = fullfile(cachedir, "topotoolbox");
else
    cachedir = options.cachedir;
end

if ~isfolder(cachedir)
    mkdir(cachedir);
end
