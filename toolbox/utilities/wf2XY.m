function [x,y] = wf2XY(wf,siz)

% Convert worldfile matrix to coordinate vectors
%
% Syntax
%
%     [x,y] = wf2XY(wf,siz)
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
%

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