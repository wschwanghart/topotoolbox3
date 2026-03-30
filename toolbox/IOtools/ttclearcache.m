function ttclearcache(options)
%TTCLEARCACHE Clear the local TopoToolbox cache directory
%
% Description
%
%     Deletes the cache directory and all the files in it. The default
%     cache directory as returned by ttcachedir will be returned unless
%     the cachedir option is provided, in which case that directory will be
%     cleared instead.
%
% Input arguments
%
%     Parameter name/value pairs
%
%     'cachedir'   char or string. If supplied, ttclearcache will use this
%     directory instead of the default one.
%
% Author: William Kearney (william.kearney[at]uni-potsdam.de)
% Date: 30 March 2026

arguments
    options.cachedir (1,1) = ttcachedir
end

if isfolder(options.cachedir)
    rmdir(options.cachedir,'s');
end


end