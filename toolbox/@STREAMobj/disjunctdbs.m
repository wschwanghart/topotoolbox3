function C = disjunctdbs(S,ixgrid)

%DISJUNCTDBS Identify mutually disjunct downstream basins in a stream network
%
% Syntax
%
%     C = disjunctdbs(S,ixgrid)
%
% Description
%
%     C = DISJUNCTDBS(S,ixgrid) partitions a set of stream pixels ixgrid
%     into groups of mutually disjunct downstream basins. The function
%     takes a stream network S (a STREAMobj) and a list of unique linear
%     indices ixgrid that lie on the stream network, and returns C, a cell
%     array in which each element contains a subset of ixgrid. Each subset
%     represents a group of pixels for which no member drains into another
%     (i.e., their basins are mutually disjunct).
%
%     The function is particularly helpful when processing nested
%     catchments defined by their outlets ixgrid. Processing subsets stored
%     in C makes sure that the drainage basins are not nested.
%
% Input arguments
%
%     S       STREAMobj
%     IXGRID  Vector of unique linear indices (into the DEM grid)
%             referring to locations on the stream network. These points
%             are typically confluences, monitoring points, or other
%             target locations.
%
% Output arguments
%
%     C       Cell array
%             C{k} contains the indices from IXGRID that form one
%             topologically disjunct group. No point in C{k} drains to
%             another point in C{k}. Together, the sets in C partition
%             IXGRID according to the downstream reachability structure
%             implied by S.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     ix = randlocs(S,100);
%     C = disjunctdbs(S,ix);
%     figure;
%     tiledlayout("flow","TileSpacing","compact");
%     for r = 1:numel(C)
%         D = drainagebasins(FD,C{r}); 
%         nexttile; 
%         imagesc(shufflelabel(D)); 
%     end
%
% See also STREAMobj, digraph, DIVIDEobj
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 21. November, 2025

arguments
    S   STREAMobj
    ixgrid {mustBeAllUnique} % must be unique
end

IX = zeros(size(S.x));
[I,b] = ismember(ixgrid,S.IXgrid);
if ~all(I)
    error('All linear indices must lie on the stream network')
end

n  = numel(ixgrid);
IX(b) = 1:n;
for r = numel(S.ix):-1:1
    if IX(S.ixc(r)) ~= 0 && IX(S.ix(r)) == 0
        IX(S.ix(r)) = IX(S.ixc(r));
    end
end

% Get a reduced stream network
I = (IX(S.ix) ~= IX(S.ixc)) & ...
    (IX(S.ixc) ~= 0);

% Transclosure (maybe use digraph/transclosure) 
A = sparse(IX(S.ix(I)),IX(S.ixc(I)),true,n,n);
A  = (speye(n)-A)\speye(n);

C  = cell(0);
ct = 0;
while nnz(A)>0
    I = sum(A,1) == 1;
    if any(I)
        ct = ct+1;
        C{ct} = ixgrid(I);
        A(I,:) = 0;
    end
end
end

% -----------------------------------------------
function mustBeAllUnique(x)
if numel(unique(x)) ~= numel(x)
    error('All values must be unique.')
end
end