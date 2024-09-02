function IX = coord2ind(FD,x,y)
%COORD2IND Convert linear indices to world coordinates
%
% Syntax
%
%     IX = coord2ind(FD,x,y)
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 31. August, 2024 


[X,Y] = wf2XY(FD.wf,FD.size);
IX = coord2ind(X,Y,x,y);