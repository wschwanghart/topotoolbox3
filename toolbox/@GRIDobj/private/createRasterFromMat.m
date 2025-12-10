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
if isequal(size(X),size(Z))
    X = X(1,:);
else
    % force row vector
    X = X(:)';
end

if isequal(size(Y),size(Z))
    Y = Y(:,1);
else
    % force column vector
    Y = Y(:);
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

% Create the referencing matrix. When created from coordinate vectors or
% matrices, the GRIDob will have a cell reference
if license('test','MAP_Toolbox')
    % Now, build a mapcellref if mapping toolbox is available.
    R = maprefcells([min(X)-cs/2 max(X)+cs/2],...
        [min(Y)-cs/2 max(Y)+cs/2],size(Z),...
        'ColumnsStartFrom','north','RowsStartFrom','west');

    % Finally check, whether
    if ~verLessThan('matlab','26')
        assert(isapprox(R.CellExtentInWorldX,csx,"tight") && ...
            isapprox(R.CellExtentInWorldY,-csy,"tight"), ...
            'Something went wrong.')
        assert(isapprox(abs(R.CellExtentInWorldX),...
                        abs(R.CellExtentInWorldY),"tight"),...
            "Cellsize in x and y direction must be the same.")
    else
        assert(isapprox_bef2024b(R.CellExtentInWorldX,csx) && ...
            isapprox_bef2024b(R.CellExtentInWorldY,-csy), ...
            'Something went wrong.')
        assert(isapprox_bef2024b(abs(R.CellExtentInWorldX),...
                        abs(R.CellExtentInWorldY)),...
            "Cellsize in x and y direction must be the same.")
    end
    wf = worldFileMatrix(R);
else
    % If mapping toolbox is not available, georef will be empty       
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


function tf = isapprox_bef2024b(a, b, opts)
%ISAPPROX   Determine approximate equality between numbers/arrays.
%
%   TF = ISAPPROX(A, B) uses default tolerances:
%       RelTol = 1e-5
%       AbsTol = 0
%
%   TF = ISAPPROX(A, B, opts) where opts is a struct with:
%       opts.RelTol  - relative tolerance
%       opts.AbsTol  - absolute tolerance

    arguments
        a
        b
        opts.RelTol (1,1) double {mustBeNonnegative} = 1e-5
        opts.AbsTol (1,1) double {mustBeNonnegative} = 0
    end

    diffVal = abs(a - b);
    tol = opts.AbsTol + opts.RelTol .* max(abs(a), abs(b));
    tf = diffVal <= tol;
end
