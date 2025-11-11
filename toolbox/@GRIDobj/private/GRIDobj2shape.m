function [MS,xy] = GRIDobj2shape(I,options)
%GRIDobj2shape Convert GRIDobj to mappolyshape or geopolyshape
%
% Syntax
%
%     [MS,xy] = GRIDobj2shape(I)
%     [MS,xy] = GRIDobj2shape(I, pn, pv, ...)
%     
% See also: GRIDobj/reclassify, GRIDobj/cropbyregion,
%           GRIDobj/GRIDobj2geotable
% 
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 7. November, 2025 

arguments (Input)
    I   GRIDobj
    options.multipart = true
    options.conn      = 8;
    options.holes     = true;
    options.isgeo     = false;
end

wf        = I.wf;
isgeo     = options.isgeo;
conn      = options.conn;

if options.holes
    holes = true;
    htxt = "holes";
else
    holes = false;
    htxt = "noholes";
end

[B,~,~,A] = bwboundaries(I.Z,conn,htxt,...
    "TraceStyle","pixeledge",...
    "CoordinateOrder","xy");

% Make sure that all polygons are closed
B = cellfun(@closePolygon,B,'UniformOutput',false);

% According to the Shapefile white paper (page 8), polygons are stored
% in a clockwise fashion except for interior polygon hole parts (i.e.
% 'donut holes'), which are stored in a counter-clockwise fashion. In
% fact, that's how polygon holes are recognized. --> ispolycw

iscw = cellfun(@(xy) ispolycw(xy(:,1),xy(:,2)),B);
% iscw returns false only, so that all boundaries are returned in cc-direction

if holes
    % This part of the code goes through all boundaries and assembles
    % enclosing and enclosed boundaries, so that boundaries of a region
    % are stored together with their holes in an nan-punctuated
    % coordinate vectors.

    G = digraph(A);
    rootnodes = indegree(G) > 0 & outdegree(G) == 0;
    D = distances(flipedge(G),find(rootnodes));

    % internal boundaries
    isinternalboundary = any(mod(D+1,2)==0,1);
    regix = find(~isinternalboundary);

    %[j,i] = find(A);
    % Identify children that are also parents and keep them in the list
    % of regions to process
    %j     = setdiff(j,i);
    %regix = setdiff((1:numel(B)),j(:));

    REG = cell(numel(regix,1));
    counter = 0;
    for k = regix
        counter = counter + 1;
        % Parents are clockwise
        if iscw(k)
            REG{counter} = [B{k}; [nan nan]];
        else
            REG{counter} = [flipud(B{k}); [nan nan]];
        end
        if (nnz(A(:,k))>0)
            for l = find(A(:,k))'
                if iscw(l)
                    REG{counter} = [REG{counter}; flipud(B{l}); [nan nan]];
                else
                    REG{counter} = [REG{counter}; B{l}; [nan nan]];
                end
            end
        end

    end
else

    REG = B;
    for k = 1:numel(REG)
        if ~iscw(k)
            REG{k} = flipud(REG{k});
            REG{k} = [REG{k};[nan nan]];
        else
            REG{k} = [REG{k};[nan nan]];
        end
    end
end

% get coordinates
XY = cellfun(@(cr)(wf*[cr-1 ones(size(cr,1),1)]')',REG,...
    'UniformOutput',false);


% if the shapes are multipart features
xy = vertcat(XY{:});
% xy(end,:) = []; % tried to remove last nan, but doesn't solve the
% problem
if isgeo
    xy = flipud(xy);
    MS = geopolyshape(xy(:,2),xy(:,1));
else
    MS = mappolyshape(xy(:,1),xy(:,2));
end

end

function xy = closePolygon(xy)

if any(xy(1,:) ~= xy(end,:))
    xy(end+1,:) = xy(1,:);
end
end