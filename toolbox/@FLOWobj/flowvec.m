function [U,V] = flowvec(FD,L)

%FLOWVEC Velocity vectors from FLOWobj
%
% Syntax
%
%     [U,V] = flowvec(FD)
%     [U,V] = flowvec(FD,L)
% 
% Description
%
%     flowvec returns the velocity components of the flow network in FD.
%     The components are normalized so that vector length is one. The
%     vector
%
% Input arguments
%
%     FD    FLOWobj
%     L     GRIDobj with length of vectors (e.g. L = flowacc(FD))
%
% Output arguments
%
%     U,V   GRIDobjs with x and y components of velocity vectors
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(fillsinks(DEM),'multi');
%     A = flowacc(FD);
%     [U,V] = flowvec(FD,sqrt(A));
%     [X,Y] = getcoordinates(DEM,'mat');
%     quiver(X,Y,U.Z,V.Z)
%
% See also: FLOWobj, FLOWobj/INFLUENCEMAP, FLOWobj/FLOWobj2M
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 31. August, 2024

arguments
    FD FLOWobj
    L  GRIDobj = GRIDobj(FD)+1
end

U = GRIDobj(FD);
[X,Y] = getcoordinates(U);
X = single(X./FD.cellsize);
Y = single(Y./FD.cellsize);
[X,Y] = meshgrid(X,Y);

DX = X(FD.ixc) - X(FD.ix);
DY = Y(FD.ixc) - Y(FD.ix);

clear X Y

if ~ismulti(FD)
    %% Single flow directions
    V = U;
    U.Z(FD.ix) = DX;
    V.Z(FD.ix) = DY;
    
else
    %% Multiple flow directions    
    U.Z   = reshape(accumarray(FD.ix,DX,[prod(U.size) 1],...
                       @mean,nan('single')),U.size);
    V     = U;
    V.Z   = reshape(accumarray(FD.ix,DY,[prod(U.size) 1],...
                       @mean,nan('single')),U.size);  
end

% Normalize vectors to have unit length
LL = sqrt(U.^2 + V.^2);
U  = U./LL;
V  = V./LL;

if nargin == 2
    U = U.*L;
    V = V.*L;
end
 