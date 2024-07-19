function [clr,xlimits] = num2rgb(x,cmap,xlimits)

%NUM2RGB Map numeric data to rgb values
%
% Syntax
%
%     [clr,xlimits] = num2rgb(x,cmap,xlimits)
%
% Description
%
%     num2rgb takes a numeric vector x and returns a rgb-triple for each
%     entry in x using the colormap cmap. The values x are mapped linearly
%     to the colors in the colormap. Values less or greater than the range
%     in xlimits are clipped to the minimum and maximum value,
%     respectively. 
%
% Input arguments
%
%     x        numeric vector
%     cmap     colormap (default = parula)
%     xlimits  minimum and maximum of the data 
%              (default = [min(x(:)) max(x(:))]). Values outside the limits
%              are clipped to the minimum or maximum value, respectively.
%
% Output arguments
%
%     clr      n*3 matrix of rgb values, where n is the number of elements
%              in x
%     xlimits  limits of x (equals the input xlimits, if provided)
%
% See also: mat2gray, gray2ind, STREAMobj/splitbyattribute
%
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 19. July, 2024


arguments
    x {mustBeNumeric}
    cmap    = parula
    xlimits (1,2) {mustBeNumeric} = [min(x(:)) max(x(:))] 
end

% Parse colormap
if isa(cmap,'char') 
    cmap = str2func(lower(cmap));
    cmap = cmap(255);
elseif isa(cmap,"function_handle")
    cmap = cmap(255);
else
    validateattributes(cmap,{'numeric'},{"ncols",3,">=",0,"<=",1},...
        'num2rgb','xlimits',3)
end


% Map x-values to colors
x   = x(:);
ind = gray2ind(mat2gray(x,xlimits),size(cmap,1));
clr = cmap(ind+1,:);
    
