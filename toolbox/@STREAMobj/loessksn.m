function [k,zhat] = loessksn(S,DEM,A,options)

%LOESSKSN Loess-smoothed river steepness
%
% Syntax
%
%     k = loessksn(S,DEM,A)
%     k = loessksn(S,DEM,A,'pn',pv,...)
%
% Description
%
%     loessksn calculates the chitransformation of the horizontal river
%     coordinate. Based on the transformed coordinate, it then calculates
%     river gradient based on a loess regression (locally estimated
%     scatterplot smoothing). The loess regression performs a locally
%     weighted linear regression and returns the regression slope as
%     estimate of ksn. The weights are derived from a tricubic weight
%     function. Note that only downstream pixels are included in the
%     estimate.
%     
% Input parameters
%
%     S       STREAMobj
%     DEM     Digital elevation model (GRIDobj or nal)
%     A       Upslope area (GRIDobj or nal) as returned by flowacc
% 
%     Parameter name/value pairs
% 
%     'a0'          Reference area {1 m^2}
%     'mn'          river concavity (theta) {0.45}
%     'ws'          window size {11}
%     'parallel'    run in parallel {false} or true
%
% Output arguments
%
%     k       node-attribute list (nal) with smoothed ksn values
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     A = flowacc(FD);
%     S = STREAMobj(FD,A>1000);
%     S = klargestconncomps(S);
%     DEM = imposemin(S,DEM);
%     k = loessksn(S,DEM,A,'parallel',false);
%     plotc(S,k)
%
% See also: STREAMobj/chitransform, STREAMobj/ksn, STREAMobj/smooth
% 
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 2. July, 2024

arguments
    S   STREAMobj
    DEM {mustBeGRIDobjOrNal(DEM,S)}
    A   {mustBeGRIDobjOrNal(A,S)}
    options.mn (1,1) {mustBeNumeric,mustBePositive} = 0.45
    options.a0 (1,1) {mustBeNumeric,mustBePositive} = 1
    options.ws (1,1) {mustBeNumeric,mustBePositive} = 11
    options.parallel (1,1) = false
end

% get node attribute list with elevation values
z = ezgetnal(S,DEM);
% get node attribute list with upstream area values
a = ezgetnal(S,A);

% ------ run in parallel -------
if options.parallel
    [CS,locb] = STREAMobj2cell(S);
    Cz = cellfun(@(ix) z(ix),locb,'UniformOutput',false);
    Ca = cellfun(@(ix) a(ix),locb,'UniformOutput',false);

    Ck = cell(numel(CS),1);
    params = options;
    params.parallel = false;
    params = namedargs2cell(params);

    parfor r = 1:numel(CS)
        Ck{r} = loessksn(CS{r},Cz{r},Ca{r},params{:});
    end

    k = getnal(S);
    for r = 1:numel(CS)
        k(locb{r}) = Ck{r};
    end
    return
end

% ------- parallel ends here -----------
    
n = numel(S.x);
M = sparse(S.ix,S.ixc,1,n,n);
M = M+speye(n);
M = M^options.ws;

d = chitransform(S,a,'a0',options.a0,'mn',options.mn);
% d = S.distance;

[i,j]  = find(M);
NEIGHS = accumarray(i,j,[n 1],@(x){x});

if nargout == 1
    k = cellfun(@(i,j) wregress(i,j),num2cell((1:n)'),NEIGHS);
else
    [k,zhat] = cellfun(@(i,j) wregress(i,j),num2cell((1:n)'),NEIGHS);
end

function [b,zhat] = wregress(ix,ixc)
    dij = d(ixc);
    if numel(dij)<=2
        b = 0;
        zhat = z(ix);
        return
    end

    wij = d(ixc)-d(ix);
    wij = wij / max(abs(wij));
    % calculate tricubic weight function
    wij = (1-abs(wij).^3).^3;
    wij = max(wij,0);
%     wij(abs(dij)>1) = 0;
    wij = wij/sum(wij);

    zij = z(ixc);
    nr  = numel(ixc);
    W   = spdiags(sqrt(wij),0,nr,nr);
    b   = (W*[ones(nr,1) dij])\(W*zij);
%     if nargout == 2
        zhat = [1 d(ix)]*b;
%     end

    b   = b(2);

end
end