function [s,dist] = stat(SW,fun,val)

%STAT Calculate statistics for swath profile object (SWATHobj)
%
% Syntax
%
%    s = stat(SW)
%    s = stat(SW,fun)
%    s = stat(SW,'prctile',val)
%    s = stat(SW,{fun1 fun2 ...})
%    [s,d] = ...
%
% Description
%
%     stat(SW) calculates the arithmetic mean across a SWATHobj for each
%     point along it.
%
%     stat(SW,fun) calculates a statistical metric across a SWATHobj for 
%     each point along it. fun specifies an aggregation function and can be
%     a handle of a function (e.g. @mean, @min, ...) that takes a vector
%     and returns a scalar, operating along the first dimension of an
%     array. Alternatively, fun can be string (e.g., 'mean', 'min', ...).
%     In this case, the text will be converted to a function handle. fun
%     can also be a 'prctile'. In this case, the function requires a third
%     input argument val, which specifies the percentile to be computed.
%     Finally, fun can be a cell array of functions. In this case, stat
%     returns a matrix s that has as many columns as there are functions.
%
% Input arguments
%
%     SW     swath object (Class: SWATHobj)
%     fun    string specifying statistical metric 
%            {'min','max','mean' (default),'range','prctile'} or anonymous
%            function (e.g. @mean). 
%
% Output
%
%     s      vector with statistic metric values along the SWATHobj
%     d      distance along swath
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
%         Wolfgang Schwanghart (schwanghart[at]uni-potsdam.de)
% Date: 11. September, 2024
    
arguments
    SW   SWATHobj
    fun  = @mean
    val  = []
end

% If multiple functions are supplied in a cell array, each will be
% calculated here
if iscell(fun)
    s = cellfun(@(f) stat(SW,f),fun,'UniformOutput',false);
    s = horzcat(s{:});
    if nargout > 1
        dist = SW.distx;
    end
    return
end

if ischar(fun) || isstring(fun)
    switch lower(fun)
        case 'prctile'
            if isempty(val)
                error('The third input argument must not be empty.')
            else
                validateattributes(val,{'numeric'},{'>',0,'<',100},'SWATHobj/stat','val',3)
            end
            calcprctile = true;
        otherwise
            calcprctile = false;
    end    
    fun = str2func(fun);
else
    calcprctile = false;
end

% Input to function
X = SW.Z;
if calcprctile
    s = fun(X,val);
else
    s = fun(X);
end
s = s';

if nargout > 1
    dist = SW.distx;
end



