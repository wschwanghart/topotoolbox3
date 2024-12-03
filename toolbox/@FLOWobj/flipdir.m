function FD = flipdir(FD,options)

%FLIPDIR Flip direction of flow
%
% Syntax
%
%     FDm = flipdir(FD)
%
% Description
%
%     flipdir reverses the direction of flow in a single or multi FLOWobj
%     so that flow moves upstream. In upstream direction, single FLOWobjs
%     become divergent and thus flipdir always returns a FLOWobj with
%     multiple flow directions.
%
%     Since formerly incoming edges in the flow network become outgoing
%     edges, the fraction along these edges needs to be computed. Moreover,
%     they need to normalized, so that they add to one for each pixel. How
%     these weights are determined is controlled by the parameter value
%     'weights'. By default, weights are determined according to the flow
%     proportions along edges calculated using flow accumulation. Hence, if
%     the pixel i has three incoming edges with 20 50 and 10. Then the
%     outgoing edges from i in the flipped direction will have values of
%     20/80, 50/80 and 10/80. 
%
% Input arguments
%
%     FD     FLOWobj
%
% Output arguments
%
%     FDm    flipped FLOWobj (multi)
%     
%     Parameter name/value pairs
%
%     'weights'  Default: []
%                If empty, then weights will be calculated according to 
%                A = flowacc(FD); w = A.Z(FD.ix).*FD.fraction
%                If weights is a GRIDobj G, then weights are calculated by 
%                w = G.Z(FD.ix).*FD.fraction
%                If weights is 'rand', then they are computed with a
%                uniformly distributed numbers
%                If weights is 'rescale', then weights are computed from
%                the normalization of the original, incoming proportions.
%                Alternatively, weights can be a vector of nonnegative that
%                has the same size as FD.ix.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(fillsinks(DEM),'multi');
%     [FD,S]  = multi2single(FD,'minarea',500);
%     FDf = flipdir(FD);
%     H = ~STREAMobj2GRIDobj(S);
%     A = flowacc(FDf,H);
%     imageschs(DEM,min(A,1000))
%
%
% See also: FLOWobj, FLOWobj/multi2single
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 3. December, 2024

arguments
    FD  FLOWobj
    options.weights =[]
end

if isempty(FD.fraction)
    FD.fraction = ones(size(FD.ix));
end

% Handle weights
if isempty(options.weights)
    A = flowacc(FD);
    w = A.Z(FD.ix).*FD.fraction;
elseif ischar(options.weights)
    switch  lower(options.weights)
        case 'rand'
           w = rand(size(FD.fraction));
        case 'rescaled'
           w = FD.fraction;
    end
elseif isa(options.weights,'GRIDobj')
    validatealignment(FD,options.weights);
    w = options.weights.Z(FD.ix).*FD.fraction;
    mustBeNonnegative(w)
elseif isnumeric(options.weights)
    if ~isequal(size(options.weights),size(FD.ix))
        error('TopoToolbox:wrongInput',...
            ['Weights must be a column vector with ' ...
              num2str(numel(FD.ix)) ' elements.'])
    end
    
    w = options.weights;
    mustBeNonnegative(w)
end

FD.type = 'multi';
ix     = FD.ix;
FD.ix  = FD.ixc;
FD.ixc = ix;

% flip ordering
FD.ix  = FD.ix(end:-1:1);
FD.ixc = FD.ixc(end:-1:1);

FD.fraction = w(end:-1:1);

FD = multi_normalize(FD);


