function L = greedylabel(D1)
%GREEDYLABEL Relabel so that adjacent regions do not share the same label
%
% Syntax
%
%     L2 = greedylabel(L1)
%
% Description
%
%     This function implements a simplified version of the map-coloring
%     algorithm. It returns a small number of labels that satisfy the
%     condition that adjacent regions do not share the same color. The
%     solution may not be optimal, though, and a smaller number of labels
%     may be theoretically possible. The function is called greedylabel
%     because it makes use of a greedy algorithm, i.e. a  problem-solving
%     approach that makes the locally optimal choice at each step with the
%     hope that these local choices will lead to a globally optimal
%     solution.
%
% Input arguments
%
%     L1     labeled grid (e.g. returned by drainagebasins)
%
% Output arguments
%
%     L2     labeled grid
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD     = FLOWobj(DEM,'single');
%     D = drainagebasins(FD);
%     L = greedylabel(D);
%     imagesc(L)
%     nlabels = max(L)
%     colormap(ttscm('davos',nlabels,[70 90]))
%
% See also: FLOWobj/drainagebasins, GRIDobj/shufflelabel
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 21. July, 2025

arguments
    D1   GRIDobj
end

% Get label values between 1 and n = the number of labels
I = D1.Z ~= 0;
flatten = @(x)x(:);
[lab,~,D1.Z(I)] = unique(flatten(D1.Z(I)));

% Sort labels so that largest regions get smallest labels
N = histcounts(flatten(D1.Z(I)),1:numel(lab)+1);
[~,ix] = sort(N,'descend');
[~,D1.Z(I)] = ismember(D1.Z(I),ix);

% Adjacency matrix of labels
D2 = dilate(D1,ones(3));
I  = D1 ~= D2;
d1 = D1.Z(I.Z);
d2 = D2.Z(I.Z);

I  = d1 == 0 | d2 == 0;
d1(I) = [];
d2(I) = [];
n  = max(D1);

% Undirected adjacency matrix
A = sparse(d1,d2,1,n,n);
A = max(A,A');

c = greedyColoring(A);
I = D1.Z ~= 0;
L = GRIDobj(D1,'single');
L.Z(I) = c(D1.Z(I));

end

function colors = greedyColoring(A)
    % Number of nodes
    n = size(A, 1);

    % Initialize colors (0 means unassigned)
    colors = zeros(1, n);

    % Assign colors to each node
    for u = 1:n
        % Find colors of adjacent nodes
        neighborColors = colors(A(u,:)>0);
        
        % Get the smallest available color not used by neighbors
        c = 1;
        while ismember(c, neighborColors)
            c = c + 1;
        end

        % Assign the color
        colors(u) = c;
    end
end
