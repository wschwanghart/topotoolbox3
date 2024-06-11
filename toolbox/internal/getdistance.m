function [cumdxy] = getdistance(x,y)

% cumulative distance along path defined by vertice coordinates
%
% Syntax
%
%      D = getdistance(x,y)
%
%
% Author: Dirk Scherler (scherler[at]@gfz-potsdam.de)
% Date: June 2013


dx = diff(x(:));
dy = diff(y(:));
dxy = hypot(dx, dy); % square root of sum of squares
cumdxy = [0; cumsum(dxy,1)];
