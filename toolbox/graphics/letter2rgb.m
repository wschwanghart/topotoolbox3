function rgb = letter2rgb(letter)

%LETTER2RGB Convert letter to rgb triple
%
% Syntax
%
%     rgb = letter2rgb(letter)
%
% Description
%
%     This function converts a letter that can be used as short-cut for color by
%     functions such as plot to a rgb triple.
%
% Example
%
%     rgb = letter2rgb('y')
%
%     returns [1 1 0]
%
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 12. November, 2024

switch lower(letter)
    case 'r'
        rgb = [1 0 0];
    case 'g'
        rgb = [0 1 0];
    case 'b'
        rgb = [0 0 1];
    case 'c'
        rgb = [0 1 1];
    case 'm'
        rgb = [1 0 1];
    case 'y'
        rgb = [1 1 0];
    case 'k'
        rgb = [0 0 0];
    case 'w'
        rgb = [1 1 1];
    otherwise
        error('Invalid color letter.')
end
        
