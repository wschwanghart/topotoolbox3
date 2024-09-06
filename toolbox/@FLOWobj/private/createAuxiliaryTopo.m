function [D,DEM,SILLS] = createAuxiliaryTopo(DEM,options)

%CREATEAUXILIARYTOPO Calculate auxiliary topography from DEM
%
% Syntax
%
%     [D,DEMF] = createAuxiliaryTopo(DEM,options)
%
% Description
%
%     Routing water and sediment flow across topography along topographic
%     gradients is straightforward as long as there are no topographic
%     sinks or flat areas. This function computes auxiliary elevations for
%     DEMs (digital elevation models) that enable to calculate flow
%     directions.
%
%
%
% See also: FLOWobj
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 6. September, 2024

arguments
    DEM   GRIDobj
    options.preprocess = 'carve'
    options.internaldrainage (1,1) = false
    options.sinks = []
    options.verbose (1,1) = false
end

% start preprocessing (e.g. carve or fill)
switch lower(options.preprocess)

    case 'fill'
        % fill will fill all internally drained basins
        % (unless specified in the sinks raster). At a
        % later stage, FLOWobj will find a route along the
        % centerline of flat areas.
        if options.verbose
            disp([char(datetime("now"))  ' -- Sink filling'])
        end

        if isempty(sinks)
            DEM = fillsinks(DEM);
        else
            DEM = fillsinks(DEM,sinks);
        end

    case 'carve'
        % carve will make use of the topography in
        % depressions to derive the most realistic flow
        % paths.
        if options.verbose
            disp([char(datetime("now"))  ' -- Sink filling'])
        end
        if isempty(options.sinks)
            DEMF = fillsinks(DEM);
        else
            DEMF = fillsinks(DEM,options.sinks);
        end

        % By default, weights are calculated as the
        % difference between the filled and the actual DEM.
        % There is also a weights option, which is a
        % GRIDobj with weights. But this is currently
        % undocumented

        D    = DEMF-DEM;
        D = D.Z;
        DEM  = DEMF;
end

% construct height graph
% After DEM filling, this code will identify flat sections,
% sills, and internally drained basins. By default
% internaldrainage is set to false. That means that FLOWobj
% will not attempt to route through the lowest region in a
% internally drained basin (e.g. a flat lake surface).
if options.internaldrainage
    [Iobj,SILLSobj,IntBasin] = identifyflats(DEM);
else
    [Iobj,SILLSobj] = identifyflats(DEM);
end

I     = Iobj.Z;
SILLS = SILLSobj.Z;

clear Iobj SILLSobj

% calculate sills for internal lake basins. These should be
% located within the lake to force convergent flows
if options.internaldrainage
    % There are various ways to get a lake center pixel
    % Here we choose the distance transform from outside
    % the lakes to the inside and take the locations as sills
    % where the distance is maximum.

    DD = bwdist(~IntBasin.Z,'e');
    STATS   = regionprops(IntBasin.Z,DD,'PixelIdxList','PixelValues');
    if ~isempty(STATS)
        for r=1:numel(STATS)
            [~,ixm] = max(STATS(r).PixelValues);
            STATS(r).MaxIntIX = STATS(r).PixelIdxList(ixm);

            I(STATS(r).PixelIdxList(1)) = false;
            SILLS(STATS(r).PixelIdxList(1)) = true;
        end

        ixm = [STATS(r).MaxIntIX];
        I(ixm) = false;
        SILLS(ixm) = true;
    end
    clear InBasin STATS

    % A slightly faster but less elegant approach
    % STATS   = regionprops(IntBasin.Z,'PixelIdxList');
    % for r = 1:numel(STATS);
    %    I(STATS(r).PixelIdxList(1)) = false;
    %    SILLS(STATS(r).PixelIdxList(1)) = true;
    %end
    % clear IntBasin STATS
end

if options.verbose
    disp([char(datetime("now")) ' -- Flat sections identified'])
end

% Some more preprocessing required. If the option carve is
% chosen, we derive here the costs to route through sinks
switch options.preprocess
    case 'carve'
        CarveMinVal = 0.1;
        % if ~isscalar(cweight)
        %     D = (D + cweight);
        %     D = linscale(D,0,100);
        % end

        % -- New version -- slightly faster but may require
        % more memory
        %                         STATS = regionprops(I,D,{'PixelIdxList','MaxIntensity','PixelValues'});
        %                         PixelValues = cellfun(@(pixval,maxval) (maxval-pixval).^tweight + CarveMinVal,...
        %                             {STATS.PixelValues},{STATS.MaxIntensity},...
        %                             'UniformOutput',false);
        %                         D(vertcat(STATS.PixelIdxList)) = cell2mat(PixelValues);

        % -- Old version
        CC = bwconncomp(I);
        tweight = 1;
        for r = 1:CC.NumObjects
            maxdepth = max(D(CC.PixelIdxList{r}));
            D(CC.PixelIdxList{r}) = (maxdepth - D(CC.PixelIdxList{r})).^tweight + CarveMinVal;
            %                             D(CC.PixelIdxList{r}) = maxdepth - D(CC.PixelIdxList{r});
            %                             D(CC.PixelIdxList{r}) = (D(CC.PixelIdxList{r})./maxdepth + CarveMinVal).^tweight;
        end
        clear CC


end

% Compute the linear index of pixels upstream of sill pixels
PreSillPixel = getPreSillPixels(DEM.Z,I,SILLS);

% Some more preprocessing if option fill is chosen. Here we
% derive the costs to route over flats
I = ~I;
switch lower(options.preprocess)
    case {'fill','none'}
        D = bwdist(I,'euclidean');
        mask = inf(FD.size,class(D));
        mask(I) = 0;
        D = (imreconstruct(D+1,mask) - D)*FD.cellsize;
    case 'carve'

end
if options.verbose
    disp([char(datetime("now"))  ' -- Weights for graydist calculated'])
end

% Here we calculate the auxiliary topography. That is, the
% cost surface seeded at socalled PreSillPixels, i.e. the
% pixel immediately upstream to sill pixels.
D(I) = inf;
D = graydist(double(D),double(PreSillPixel),'q') + 1;
D(I) = -inf;

if options.verbose
    disp([char(datetime("now")) ' -- Auxiliary topography in flats calculated'])
end
