function [mn,results] = mnoptim(S,DEM,A,options)

%MNOPTIM Bayesian optimization of the mn ratio
%
% Syntax
%   
%     [mn,results] = mnoptim(S,DEM,A)
%     [mn,results] = mnoptim(S,z,a)
%     [mn,results] = mnoptim(...,pn,pv,...)
%
% Description
%
%     mnoptim uses Bayesian optimization to find an optimal mn-ratio by
%     minimizing a cross-validation loss. Cross-validation is performed
%     using repeated 2-fold validation on individual stream networks
%
% Input arguments
%
%     S      STREAMobj
%     DEM    Digital elevation model (GRIDobj)
%     A      upstream area (in pixels) as returned by the function flowacc
%            (GRIDobj)
%     z      node-attribute list of elevations
%     a      node-attribute list of upstream areas (nr of upstream pixels)
%
% Parameter-name value pairs
%
%     'optvar'    {'mn'},'m','n','m&n'
%                 The variable to be optimized. 'mn' optimizes the mn ratio,
%                 'm' finds an optimal value of m for a given value of n.
%                 'n' finds an optimal value of n for a given value of m.
%                 'm&n' find an optimal value of m and n.
%     'n'         value of n, only applicable if 'optvar' is 'm'. Default
%                 value is 1.
%     'm'         value of m, only applicable if 'optvar' is 'n'. Default
%                 value is 0.5. 
%     'nrange'    search range of 'n'. Only applicable if 'optvar' is 'n' or
%                 'm&n'. Default is [0.8 1.5].
%     'mrange'    search range of 'm'. Only applicable if 'optvar' is 'm'.
%                 'm&n'.Default is [0.3 1].
%     'mnrange'   search range of 'mn'. Only applicable if 'optvar' is 'm'.
%                 Default is [.1 1].
%     'a0'        reference area (see chitransform). Default is 1e6 m^2.
%     'lossfun'   function to be minimized. Default is 'rmse'. Others are
%                 'linear' and 'coeffdeterm'.
%     'UseParallel' {true} or false.
%     'crossval'  {true} or false.   
%
% Other bayesopt parameter name/value pairs (see help bayesopt for details)
%
%     'MaxObjectiveEvaluations'
%
%
% Output arguments
%
%     mn       table with results (best point)
%     results  BayesianOptimization object. Continue with optimization
%              with the command resultsnew = resume(results).
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     A  = flowacc(FD);
%     S  = STREAMobj(FD,'minarea',1e6,'unit','map');
%     S  = removeedgeeffects(S,FD,DEM);
%     [mn,results] = mnoptim(S,DEM,A,'optvar','mn','crossval',true);
%     
% See also: chiplot, chitransform, slopearea
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 17. June, 2024

arguments
    S     STREAMobj
    DEM {mustBeGRIDobjOrNal(DEM,S)}
    A   {mustBeGRIDobjOrNal(A,S)}
    options.mnrange (1,2) {mustBeNumeric,mustBePositive} = [.1 1]
    options.optvar {mustBeTextScalar} = 'mn'
    options.n (1,1) {mustBeNumeric,mustBePositive} = 1
    options.m (1,1) {mustBeNumeric,mustBePositive} = 0.5
    options.nrange (1,2) {mustBeNumeric,mustBePositive} = [0.8 1.5]
    options.mrange (1,2) {mustBeNumeric,mustBePositive} = [0.3 1.2]
    options.a0 (1,1) {mustBeNumeric,mustBePositive} = 1
    options.lossfun = 'rmse'
    options.plot  (1,1) = false
    options.crossval (1,1) = true
    options.verbose (1,1) = 0
    options.UseParallel (1,1) = true
    options.MaxObjectiveEvaluations (1,1) {mustBeInteger,mustBePositive} = 30
end


a0 = options.a0;

% get node attribute list with elevation values
z = ezgetnal(S,DEM);

% get node attribute list with flow accumulation values
a = ezgetnal(S,A);
clear A DEM

% set the base level of all streams to zero
zbase = z;
for r = numel(S.ixc):-1:1
    zbase(S.ix(r)) = zbase(S.ixc(r));
end
z = z-zbase;


%% -----------------------------------------------------------------------

if options.crossval
    [~,locb] = STREAMobj2cell(S);
    cv = true;
    % Number of connected components
    nrcc = numel(locb);
    if nrcc == 1
        cv = false;
        error('TopoToolbox:mnoptim',['Cross-validation not possible. There is only \n' ...
                                       'one connected component in the stream network.\n' ...
                                       'Set option ''crossval'' to false.']);
    end
else
    cv = false;
end

% get lossfunction
if ischar(options.lossfun)
    switch options.lossfun
        case 'rmse'
            lossfun = @(x,xhat)mean(sum((x-xhat).^2));
        case 'linear'
            lossfun = @(x,xhat)var((x+1)./(xhat+1));
        case 'linearcc'
            cc = conncomps(S);
            lossfun = @(x,xhat) sum(accumarray(cc,(x+1)./(xhat+1),[],@std)).^2;
%             lossfun = @(x,xhat)var((x+1)./(xhat+1));
            
        case 'coeffdeterm'
            take2elem = @(x) x(2);
            lossfun = @(x,xhat)1-(take2elem(corrcoef(x,xhat))).^2;
        otherwise
            error('unknown loss function')
    end
else
    lossfun = options.lossfun;
end

% Bayes
optvar = options.optvar;
isdeterm = ~cv;
pp   = gcp('nocreate');
opts = {'IsObjectiveDeterministic', isdeterm, ...
        'Verbose', options.verbose,...
        'UseParallel',~(isempty(pp) & ~options.UseParallel),...
        'MaxObjectiveEvaluations',options.MaxObjectiveEvaluations,...
        };
        
switch optvar
    case 'mn'
        mn      = optimizableVariable('mn',options.mnrange);    
        results = bayesopt(@(x) fun(x),mn,opts{:});
    case 'm'
        m       = optimizableVariable('m',options.mrange);
        n       = options.n;
        results = bayesopt(@(x) fun(x),m,opts{:});    
    case 'n'
        n       = optimizableVariable('n',options.nrange);
        m       = options.m;
        results = bayesopt(@(x) fun(x),n,opts{:});
    case 'm&n'
        m       = optimizableVariable('m',options.mrange);
        n       = optimizableVariable('n',options.nrange);
        results = bayesopt(@(x) fun(x),[m n],opts{:});
               
end

mn      = bestPoint(results);



% loss function
function lss = fun(x)


switch optvar
    case 'mn'
        c = chitransform(S,a,'mn',x.mn,'a0',a0);    
    case 'm'
        c = chitransform(S,a,'mn',x.m/n,'a0',a0);  
    case 'n'
        c = chitransform(S,a,'mn',m/x.n,'a0',a0);  
    case 'm&n'
        c = chitransform(S,a,'mn',x.m/x.n,'a0',a0);  
end

if ~cv    
    % no cross-validation
    b = z\c;
    zhat = b*c;
    lss = lossfun(z/max(z),zhat/max(zhat));
else
    % cross validation
    LSS = zeros(5,1);
    for iter = 1:numel(LSS)
        conncomptrain = randperm(nrcc,ceil(nrcc/2));
        CC            = false(nrcc,1);
        CC(conncomptrain) = true;
        ixtrain       = vertcat(locb{CC});
        b = z(ixtrain)\c(ixtrain);
        ixval   = vertcat(locb{~CC});
        zhat = b*c(ixval);
        zt   = z(ixval);
        LSS(iter)  = lossfun(zt/max(zt),zhat/max(zhat));
%         LSS(iter)  = lossfun(zt/prctile(zt,75),zhat/prctile(zhat,75));
    end
    lss = mean(LSS);
end
    

end
end

