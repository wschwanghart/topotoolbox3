function clr = ttclr(feature)

%TTCLR Colors for common map features
%
% Syntax
%
%     clr = ttclr(feature)
%
% Description
%
%     ttclr returns a number of rgb values for common usage. 
%
% Input arguments
%
%     feature     'lake','lakeoutline','river','glacier','desert', or 
%                 'meadow'
%
% Output arguments
%
%     clr         rgb-triple 
%
% See also: ttscm, ttcmap
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 12. October, 2024

arguments
    feature {mustBeMember(feature,{'lake','lakeoutline',...
                                   'river','glacier',...
                                   'desert','meadow'})}
end

switch lower(feature)
    case 'lake'
        clr = [165, 191, 221]/255;
    case {'lakeoutline','river'}
        clr = [0 0.4470 0.7410];
    case 'glacier'
        % https://support.esri.com/en/technical-article/000010027
        % 10% Gray
        clr = [225 225 225]/255;
    case 'desert'
        % https://support.esri.com/en/technical-article/000010027
        % Sahara sand
        clr = [255 235 190]/255;
    case 'meadow'
        clr = [48 186 143]/255;
end