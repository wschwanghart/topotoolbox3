function [v,d,xy] = applyfun(SW,fun,A)

%APPLYFUN Apply function to extract values along swath
%
% Syntax
%
%     [v,d] = applyfun(SW,fun)
%     [v,d] = applyfun(SW,fun,A)
%
% Description
%
%     This functions applies a function fun to summarize/aggregate values
%     along the swath SW. fun must be a function that takes a 2D-matrix and
%     computes values along the first dimension of an array (e.g. @mean).
%     fun can also be a cell array of functions.
%
% Input arguments
%
%     SW    SWATHobj
%     fun   function or character or cell array of function
%     A     GRIDobj
%
% Output arguments
%
%     v     values mapped to trace
%     d     distance along trace
%     xy    trace coordinates
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     x = [382209 389559]';
%     y = [3790773 3799443]';
%     SW = SWATHobj(DEM,x,y,'width',1000);
%     [v,d] = applyfun(SW,@(x) quantile(x,[.1:.1:.9]));
%     plot(d,v)
%
% See also: SWATHobj, SWATHobj/plotdz
%           
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 27. March, 2026

arguments
    SW   SWATHobj
    fun  
    A = []
end

if ~isempty(A)
    SW = mapswath(SW,A);
end

% Place function in cell array
if ~iscell(fun)
    fun = {fun};
end

% Iterate through cell array of functions
for r = 1:numel(fun)
    % Convert string or character array to anonymous function
    if ischar(fun{r}) || isstring(fun{r})
        f = str2func(fun{r});
    else
        f = fun{r};
    end

    if r == 1
        v = f(SW.Z);
    else
        v = vertcat(v,f(SW.Z)); %#ok<AGROW>
    end
end

v = v';

if size(v,1) ~= size(SW.Z,2)
    error('TopoToolbox:error',...
        'Function output inconsistent with size of trace.')
end
d  = SW.distx;
xy = SW.xy;


