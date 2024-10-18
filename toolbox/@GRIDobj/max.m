function [m,varargout] = max(varargin)

%MAX Maximum value in GRIDobj
%
% Syntax
%
%     m = max(DEM)
%     [m,ix] = max(DEM)
%     DEMm   = max(DEM,DEM2)
%
% Description
%
%     Returns the maximum value in a GRIDobj
%
% Input arguments
%
%     DEM   GRIDobj
%     DEM2  GRIDobj (must spatially align with DEM)
%
% Output arguments
%
%     m     maximum value
%     ix    location of maximum value (linear index in GRIDobj)
%     DEMm  GRIDobj with elementwise minimas
%
% See also: GRIDobj/min
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 18. October, 2024


funname = 'max';

narginchk(1, 2)
if nargin == 1
    nargoutchk(0, 2)
    [m,varargout{1}] = max(varargin{1}.Z(:));
    
elseif nargin == 2
    nargoutchk(0, 1)
    
    % which one of the two input arguments is a GRIDobj
    isGRIDobj = cellfun(@(x) isa(x,'GRIDobj'),varargin);
    
    % if both
    if all(isGRIDobj)
        validatealignment(varargin{1},varargin{2}); 
        M = builtin(funname,varargin{1}.Z,varargin{2}.Z);
    else
        varargin = varargin(~isGRIDobj + 1);
        if ~isscalar(varargin{2})
            validatealignment(varargin{1},varargin{2}); 
        end
        M = builtin(funname,varargin{1}.Z,varargin{2});
    end

    m = varargin{1};
    m.name = funname;
    m.Z = M; 
    m.zunit = '';
end
