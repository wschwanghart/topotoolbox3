function D = meanflowdist(FD,S,options)

%MEANFLOWDIST Average flow distance based on multiple flow directions
%
% Syntax
%
%     D = meanflowdist(FD)
%     D = meanflowdist(FD,S)
%     D = meanflowdist(__,'direction','downstream' or 'upstream')
%
% Description
%
%     This function calculates the average distance (average pathlength)
%     that is required to reach each pixel in a DEM. The function can be
%     used with multiple and single flow directions. If a STREAMobj is
%     supplied (e.g. as output of the function multi2single), then average
%     downstream distances will be calculated for pixels upstream of the
%     flow network. Values in D at stream locations then refer to the
%     average length of the hillslope upstream to the channel location
%     (excluding upstream areas along the stream network).
%
% Input arguments
%
%     FD    FLOWobj (multi or single)
%     S     STREAMobj.
%
%     Parameter name/value pairs
%
%     'direction'   {'downstream'} or 'upstream'
%     'seed'        GRIDobj or linear index for seed locations from which 
%                   distances will be calculated 
%
% Output arguments
%
%     D     GRIDobj with average donwstream distances
%
% Example 1
%
%     [X,Y] = meshgrid(-5:0.01:5);
%     % Bivariate Gaussian distribution 
%     Z = 1/(2*pi) * exp(-(1/2)*(X.^2 + Y.^2));
%     DEM = GRIDobj(X,Y,Z);
%     FD  = FLOWobj(DEM,'dinf');
%     D   = meanflowdist(FD);
%     subplot(1,2,1)
%     contour(D)
%     axis image
%     subplot(1,2,2)
%     plot(D.Z(:),DEM.Z(:),'.');
%     xlabel('Mean flow distance'); ylabel('Elevation')
%
% Example 2
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'multi');
%     [FD,S] = multi2single(FD,'minarea',1000);
%     D = meanflowdist(FD,S);
%     imagesc(D)
%
% Example 3
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'multi');
%     D = meanflowdist(FD,'seed',663859);
%     imageschs(DEM,D)
%
% Example 4
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'multi');
%     [x,y] = randomsample(DEM,100);
%     ix = coord2ind(DEM,x,y);
%     D = meanflowdist(FD,'seed',ix);
%     imageschs(DEM,D)     
%
% See also: FLOWOBJ, FLOWobj/flowdistance, FLOWobj/flowconvergence,
%           FLOWobj/fHBR, FLOWobj/fHAND, FLOWobj/dbentropy
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 17. March, 2026

arguments
    FD   FLOWobj
    S    = []
    options.A {mustBeGRIDobjOrEmpty} = []
    options.direction {mustBeMember(options.direction,{'upstream',...
        'downstream'})} = 'downstream'
    options.seed = [];
end


switch options.direction
    case 'downstream'
        ix  = FD.ix;
        ixc = FD.ixc;
        p   = FD.fraction;

        if isempty(p)
            % If FD stores single flow directions, there is no fraction
            % stored. Thus, we generate a vector of ones.
            p = ones(size(ix));
        end

        % if a seed is supplied, we will restrict the flow network to the
        % subgraph sourced at the seed pixel
        if ~isempty(options.seed)
            I = influencemap(FD,options.seed);
            I = I.Z(ix);
            ix = ix(I);
            ixc = ixc(I);
            p   = p(I);
            clear I
        end
    
        
        if isempty(options.A)
            A = flowacc(FD);
        else
            A = options.A;
            validatealignment(FD,A);
        end

        % Create a vector of edge weights (fraction*upslopearea)
        p   = p.*A.Z(ix);

        % If a STREAMobj is supplied, then node-to-node links between channel
        % locations are removed from the FLOWobj edge-list
        if ~isempty(S)
            I = ~ismember([ix ixc],[S.IXgrid(S.ix) S.IXgrid(S.ixc)],"rows");
            ix = ix(I);
            ixc = ixc(I);
            p   = p(I);
        end

        % Calculate node-to-node distance
        [X,Y] = getcoordinates(A,'matrix');
        d     = hypot(X(ix)-X(ixc),Y(ix)-Y(ixc));

        % Preallocate distance
        D   = zeros(FD.size);

        % FR sums the fractions of incoming edges of each cell. These are required
        % to normalize the weighted averages at incoming nodes.
        FR = zeros(FD.size);
        for r = 1:numel(ix)
            FR(ixc(r)) = FR(ixc(r)) + p(r);
        end

        % The average distance of a pixel (ixc) to its upstream neighbors is calculated
        % from the sum of distances of the upstream neighbors plus their edge distance
        % weighted by the area-fraction. To obtain the average, the summed value needs to
        % be normalized using the sum of weights of the incoming edges.
        %
        for r = 1:numel(ix)
            D(ixc(r)) = p(r).*(D(ix(r)) + d(r))/FR(ixc(r)) + D(ixc(r));
        end
        % Convert to a GRIDobj
        D = GRIDobj(FD,D);

    case 'upstream'
        % In upstream direction, we just need to flip the directions in FD
        FD = flipdir(FD,"weights","rescaled");
        % In case, a STREAMobj is supplied, we will calculate the distances
        % from the stream network. This means, we will simply delete edges
        % in the flow network that belong to the river network.
        if ~isempty(S)
            I = ~ismember(FD.ix,S.IXgrid);
            FD.ix = FD.ix(I);
            FD.ixc = FD.ixc(I);
            FD.fraction = FD.fraction(I);
        end

        % Now we simply need to call meanflowdist in downstream direction
        D = meanflowdist(FD,'direction','downstream');
end

