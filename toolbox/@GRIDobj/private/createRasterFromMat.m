function [Z,R,wf] = createRasterFromMat(X,Y,Z)

%CREATERASTERFROMMAT Creates a grid from coordinate and value arrays
%
% Syntax
%
%     [Z,R,wf] = createRasterFromMat(X,Y,Z)
%
%

arguments
    X {mustBeNumeric}
    Y {mustBeNumeric}
    Z {mustBeEqualSize(X,Y,Z),mustBeNumericOrLogical}
end


% Create coordinate vectors if matrices are supplied
if ismatrix(X)
    X = X(1,:);
end

if ismatrix(Y)
    Y = Y(:,1);
end

% The convention in a GRIDobj is that columns start from north, and rows
% start from west.
% Check X first
if X(2) > X(1)
    % Rows start from west.
else
    % Rows start from east -> need to flip coordinates and Z-array
    X = X(end:-1:1);
    Z = Z(:,end:-1:1);
end

if Y(1) > Y(2)
    % Columns start from north.
else
    % Columns start from west -> need to flip coordinates and Z-array
    Y = Y(end:-1:1);
    Z = Z(end:-1:1,:);
end

% Now check, whether coordinate vectors are uniformly spaced
[tfx,csx] = isuniform(X);
[tfy,csy] = isuniform(Y);

cs = abs(csx);

assert(tfx && tfy, 'Coordinate vectors/matrices must be uniformly spaced.')
%assert(csx == -csy, 'The step size (cellsize) along x and y must be the same.')

% Create the referencing matrix
if license('test','MAP_Toolbox')
    % Now, build a mapcellref if mapping toolbox is available.
    R = maprefcells([min(X)-cs/2 max(X)+cs/2],...
        [min(Y)-cs/2 max(Y)+cs/2],size(Z),...
        'ColumnsStartFrom','north','RowsStartFrom','west');

    % Finally check, whether
    assert(R.CellExtentInWorldX == csx && R.CellExtentInWorldY == -csy, ...
        'Something went wrong.')
    wf = worldFileMatrix(R);
else
    R = [];
    wf = [csx    0    min(X);
        0      csy  max(Y)];

end


end

function mustBeEqualSize(x,y,z)
% Checks whether the size of x, y and z are compatible

eid = 'TopoToolbox:Size:notEqual';
% x and y are vectors
if isvector(x) && isvector(y)
    if ~isequal(size(z,1), numel(y))
        msg = 'The length of the vector Y must equal the number of rows in Z.';
    end
    if ~isequal(size(z,2), numel(x))
        msg = 'The length of the vector X must equal the number of columns in Z.';
    end
else
    if ~isequal(size(x),size(z)) || ~isequal(size(y),size(z))
        msg = 'The size of the matrices X, Y and Z must be equal.';
    end
end

if exist('msg','var')
    error(eid,msg)
end

end