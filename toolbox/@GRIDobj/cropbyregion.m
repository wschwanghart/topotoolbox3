function [C,BL] = cropbyregion(L,I,CL)

%CROPBYREGION Extract and crop based on a GRIDobj with labelled regions
%
% Syntax
%
%     C = cropbyregion(L)
%     [C,BL] = cropbyregion(L)
%     [CDEM,BL] = cropbyregion(L,DEM)
%     [CDEM,BL] = cropbyregion(L,DEM,doclip)
%     CDEM = cropbyregion(DEM,BL)
%     CDEM = cropbyregion(DEM,BL,C)
%
% Description
%
%     This function extracts and crops a cell array of GRIDobj C based on a
%     GRIDobj with labelled regions (e.g. computed with the function
%     drainagebasins). Calling cropregion with a single input returns the
%     cell array C which is a n-by-1 cell array with logical GRIDobjs where
%     n is the number of regions in L. The second output BL is a n-by-1
%     cell array. Each element contains the subscripts of the block limits
%     [top bottom left right].
%
%     If a GRIDobj DEM is provided as second input, cropbyregion will
%     extract and crop the DEM based on the regions in L. clip is a logical
%     scalar. If true (default), values outside the regions in L will be
%     set to nan.
%
%     BL can be used as second argument to extract a cell array of cropped 
%     GRIDobjs from the GRIDobj DEM. If the argument C is included as third
%     argument, the DEM will be clipped to the regions in C. 
%
%     The function requires that L has a an associated georeferencing
%     object (L.georef). GRIDobjs in C will be expanded if regions consist
%     of single pixels. The minimum size is two rows and two columns.
%
% Input arguments
%
%     L       GRIDobj with labelled regions
%     DEM     GRIDobj to be cropped based on regions in L
%     doclip  {true} or false. If true, the DEMs in CDEM will be cropped
%             and clipped
%     BL      Cell array with block limits (used for cropping)
%     C       Cell array with cropped regions (used for clipping)
%
% Output arguments
%
%     C    cell array
%     BL   Block limits
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S   = STREAMobj(FD,'minarea',10000);
%     D   = drainagebasins(FD,S);
%     [C,BL] = cropbyregion(D);
%     CDEM = cropbyregion(DEM,BL,C);
%     tiledlayout
%     for r = 1:numel(CDEM); 
%           nexttile; 
%           imageschs(CDEM{r}); 
%     end
%
% See also: GRIDobj/crop, GRIDobj/clip, regionprops
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 6. November, 2025

arguments
    L    GRIDobj
    I   = []
    CL   = true
end

if isempty(I) || isa(I,"GRIDobj")

    if ~isempty(I)
        validatealignment(L,I)
    end
    validateattributes(CL,["numeric","logical"],"scalar","GRIDobj\cropbyregion","clip",3)
    if CL & ~isempty(I)
        if isfloat(I.Z)
            clipval = nan(1,class(I.Z));
        else
            clipval = zero(1,class(I.Z));
        end
    end

    % Obtain images and bounding boxes of regions
    stats = regionprops(L.Z,{'Area','BoundingBox','Image'});

    % Size of the grid
    siz   = L.size;

    % Compute bounding box so that it refers to pixel subscripts
    BB    = {stats.BoundingBox};
    BB    = cellfun(@(bb) [floor(bb(1:2))+1 max(bb(3:4)-1,0)],BB,...
        'UniformOutput',false);

    % Referencing object (later required by cropToBlock)
    R   = L.georef;

    % Preallocate
    C = cell(numel(stats),1);
    BL = cell(numel(stats),1);

    for r=1:numel(stats)
        if stats(r).Area < 1
            C{r} = [];
            BL{r} = [];
        else
            colLimits = [BB{r}(1)  BB{r}(1)+BB{r}(3)];
            rowLimits = [BB{r}(2)  BB{r}(2)+BB{r}(4)];

            % If any of the limits is only one pixel wide, the limits need to
            % be expanded.
            if colLimits(1) == colLimits(2)
                nrows = size(stats(r).Image,1);
                if colLimits(1) < siz(2)
                    colLimits(2) = colLimits(2)+1;
                    % add an extra column to the right
                    stats(r).Image = [stats(r).Image false(nrows,1)];
                else
                    colLimits(1) = colLimits(1)-1;
                    % add an extra column to the left
                    stats(r).Image = [false(nrows,1) stats(r).Image];
                end
            end
            if rowLimits(1) == rowLimits(2)
                ncols = size(stats(r).Image,2);
                if rowLimits(1) < siz(1)
                    rowLimits(2) = rowLimits(2)+1;
                    % add an extra row at the bottom
                    stats(r).Image = [stats(r).Image; false(1,ncols)];
                else
                    rowLimits(1) = rowLimits(1)-1;
                    % add an extra row at the top
                    stats(r).Image = [false(1,ncols); stats(r).Image];
                end
            end

            RC    = cropToBlock(R,rowLimits,colLimits);
            if nargin == 1
                C{r}  = GRIDobj(stats(r).Image,RC);
            else
                im    = I.Z(rowLimits(1):rowLimits(2),colLimits(1):colLimits(2));
                C{r}  = GRIDobj(im,RC);
                if CL
                    C{r}.Z(~stats(r).Image) = clipval;
                end
            end
            BL{r} = [rowLimits colLimits];

        end
    end

else
    BL = I;
    % Blocklimits
    validateattributes(BL,"cell",{},"GRIDobj\cropbyregion","BL",2)
    % Clip rasters
    validateattributes(CL,"cell",{},"GRIDobj\cropbyregion","C",3)

    C = cell(numel(BL),1);
    for r = 1:numel(BL)
        IX = sub2ind(L.size,BL{r}(1:2),BL{r}(3:4));
        C{r} = crop(L,IX);

        if iscell(CL)
            C{r} = clip(C{r},CL{r});
        end
    end
end

