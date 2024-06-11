function [x,y] = wf2XY(wf,siz)

%WF2XY Convert worldfile matrix to coordinate vectors
%
% Syntax
%
%     [x,y] = wf2XY(wf,siz)
%
% Description
%
%     wf2XY creates coordinate vectors from a worldFileMatrix wf and the
%     size of the array. 
%
% Input arguments
%
%     wf    worldfile matrix
%     siz   size of the grid
%
% Output arguments
%
%     x     x-coordinates
%     y     y-coordinates
%
% Example
%
%     wf = [30   0   0.37632865e6; ...
%            0 -30   3.80790282e6]
%     siz = [643 1197];
%     [x,y] = wf2XY(wf,siz);
%
% See also: GRIDobj/coord2ind, GRIDobj/getcoordinates
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 30. May, 2024

arguments
    wf {mustBeWorldFile}
    siz {mustBeNumeric,mustBePositive}
end

nrrows = siz(1);
nrcols = siz(2);

if wf(2) == 0 % if grid is rectilinear
    x = wf(1,1)*(0:nrcols-1)' + wf(1,3);
    y = wf(2,2)*(0:nrrows-1)' + wf(2,3);
else % with rotation
    % intrinsic coordinate matrices
    [Xi,Yi] = meshgrid(0:nrcols-1,0:nrrows-1);
    T = wf(1:2,1:2)*[Xi(:)'; Yi(:)'];
    x = reshape(T(1,:)',siz);
    y = reshape(T(2,:)',siz);

    x = x+ wf(1,3);
    y = y+ wf(2,3);
end
end

function mustBeWorldFile(x)

if ~isequal(size(x),[2 3])
    error('TopoToolbox:WrongInput','WorldFileMatrix must have size [2 3]');
end
end