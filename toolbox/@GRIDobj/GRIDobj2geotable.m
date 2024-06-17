function [GT,x,y] = GRIDobj2geotable(DB,options)

%GRIDOBJ2GEOTABLE Convert categorical GRIDobj to geotable (polygon) 
%
% Syntax
%
%     GT = GRIDobj2geotable(DB)
%
% Description
%
%     GRIDobj2geotable takes a GRIDobj with categorical values and stores
%     the outlines of each region with the same values as geotable. 
%
% Input arguments
%
%     DB    GRIDobj with categorical values (this is not meant that the
%           underlying type is categorical but that there are contiguous
%           regions with the same values. These need not to be integer
%           values.
% 
%     Parameter name/value pairs
%
%     'simplify'     {false} or true (NOT WORKING YET)
%     'tol'          0 (NOT WORKiNG YET)
%     'excludezero'  {true} or false. If true, regions with zero value will
%                    not be polygonized.
%     'conn'         8 or 4 connectivity. 8 is default.
%     'holes'        {true} or false
%     'parallel'     {true} or false. If true, function is run in parallel
%     'geographic'   {false} or true. If true, the returned geotable will
%                    store geographic coordinates (lat-lon).
%
% Output arguments
%
%     GT      geotable
%     x,y     nan-punctuated coordinate vectors of the polygon outlines
%     
% Example
%
%     DEM = readexample('taalvolcano');
%     I = identifyflats(DEM);
%     I.Z = bwareaopen(I.Z,20);
%     DEM = clip(DEM,~I);
%     imageschs(DEM)
%    
%     C = reclassify(DEM,'equalint',5);
%     GT = GRIDobj2geotable(C);
%     geoplot(GT)
%
% See also: GRIDobj/reclassify
% 
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 17. June, 2024 

arguments (Input)
    DB   GRIDobj
    options.simplify = false
    options.tol      = 0
    options.multipart = true
    options.excludezero = true
    options.conn      = 8;
    options.holes    = true;
    options.parallel = true;
    options.geographic = false;
end

EXCL.Z = isnan(DB.Z);
if options.excludezero
    EXCL.Z = EXCL.Z | DB.Z == 0;
end


[uniqueval,~,DB.Z(~EXCL.Z)] = unique(DB.Z(~EXCL.Z));
nuniqueval = numel(uniqueval);

wf        = DB.wf;
conn      = options.conn;
multipart = options.multipart;
prj       = DB.georef.ProjectedCRS;
geographic = options.geographic;

if options.holes
    holes = true;
    htxt = "holes";
else
    holes = false;
    htxt = "noholes";
end

if nargout > 1
    XYC = cell(nuniqueval,1);
    docoords = true;
else
    docoords = false;
end

CMS = cell(nuniqueval,1);

% Run in parallel?
if ~options.parallel
    poolsize = 0;
else
    p = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(p)
        poolsize = 0;
    else
        poolsize = p.NumWorkers;
    end
end

%parfor (uv = 1:nuniqueval,poolsize)
for uv = 1:nuniqueval

    I = DB == uniqueval(uv);
    val = uniqueval(uv);

    [B,~,~,A] = bwboundaries(I.Z,conn,htxt,"TraceStyle","pixeledge",...
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

    if geographic
        [X,Y] = cellfun(@(xy) projinv(prj,xy(:,1),xy(:,2)),XY,...
            'UniformOutput',false);
        XY = cellfun(@(x,y) [x y],X,Y,'UniformOutput',false);
    end

    if docoords
        XYC{uv} = vertcat(XY{:});
    end

    % create polyshapes
    if multipart
        % if the shapes are multipart features
        xy = vertcat(XY{:});
        % xy(end,:) = []; % tried to remove last nan, but doesn't solve the
                          % problem
        if geographic
            MS = geopolyshape(xy(:,1),xy(:,2));
        else
            MS = mappolyshape(xy(:,1),xy(:,2));
            MS.ProjectedCRS = prj;
        end
        MS = table(MS,val,'VariableNames',{'Shape','Value'});
    else
        % else if shapes are single part features
        MS = [];
        for r = 1:numel(XY)
            xy = XY{r};
            % xy(end,:) = []; % tried to remove last nan, but doesn't solve the
                              % problem
            if geographic
                S = geopolyshape(xy(:,1),xy(:,2));
            else
                S = mappolyshape(xy(:,1),xy(:,2));
                S.ProjectedCRS = prj;
            end

            if r == 1
                MS = table(S,1,r,'VariableNames',{'Shape','Value','Part'});
            else
                MS = [MS; table(S,1,r,'VariableNames',{'Shape','Value','Part'})];
            end
        end
    end

    CMS{uv} = MS;
end

GT = vertcat(CMS{:});
if nargout > 1
    XYC = vertcat(XYC{:});
    x = XYC(:,1);
    y = XYC(:,2);
end
end


function xy = closePolygon(xy)

if xy(1,:) ~= xy(end,:)
    xy(end+1,:) = xy(1,:);
end
end