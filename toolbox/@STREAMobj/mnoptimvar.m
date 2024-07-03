function [mn,cm,zm,zsd] = mnoptimvar(S,DEM,A,options)

%MNOPTIMVAR Optimize the mn ratio using minimum variance method
%
% Syntax
%
%     mn = mnoptimvar(S,DEM,A)
%     mn = mnoptimvar(S,z,a)
%     mn = mnoptimvar(...,pn,pv,...)
%     [mn,cm,zm,zsd] = ...
%
% Description
%
%     mnoptimvar finds an optimal value of the mn-ratio by minimizing the
%     variance of elevation values conditional on chi. The procedure places
%     nodes of the river network into chi-distance bins and calculates the
%     variance of elevation values in each bin. The objective function is
%     the weighted mean of these variances.
%
% Input arguments
%
%     S     STREAMobj
%     DEM   Digital elevation model
%     A     Flow accumulation grid (as returned by flowacc)
%     z     node-attribute list of elevations
%     a     node-attribute list of flow accumulation
%
%     Parameter name/value pairs
%
%     'varfun'   function to calculate elevation variability in chi-distance 
%                bins. Default is the interquartile range (@iqr), but any 
%                other measure of variability may be suitable, e.g. @var,
%                @std, @mad, @robustcov, @range.
%     'mn0'      starting value for mn {0.5}.
%     'distbins' nr of distance bins {100}.
%     'minstream' minimum number of streams that should be included in
%                calculation a variance value {2}.
%     'funoptim' {'fminsearch'} (so far only choice)
%     'optimize' {true} or false. If false, then the function won't perform
%                an optimization, but will use mn0 as mn-ratio.
%     'plot'     {true} or false
%     'zerobaselevel' {false} or true. Set all outlet elevations to zero.
%     'a0'       reference area {1e6}
%     
% Output arguments
%
%     mn    mn-ratio
%     cm    binned chi values
%     zm    average elevation in chi bins
%     zsd   standard deviation of elevations in chi bins
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     A        = flowacc(FD);
%     S   = STREAMobj(FD,'minarea',1000);
%     S   = klargestconncomps(S);
%     mn  = mnoptimvar(S,DEM,A);
%     
%
% See also: STREAMobj, STREAMobj/mnoptim
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 17. June, 2024

arguments
    S     STREAMobj
    DEM {mustBeGRIDobjOrNal(DEM,S)}
    A   {mustBeGRIDobjOrNal(A,S)}
    options.varfun = @iqr
    options.distbins (1,1) {mustBePositive,mustBeInteger} = 100
    options.minstreams (1,1) {mustBePositive,mustBeInteger} = 2
    options.mn0 (1,1) {mustBePositive,mustBeNumeric} = 0.5
    options.funoptim = 'fminsearch'
    options.plot (1,1) = true
    options.optimize (1,1) = true
    options.zerobaselevel (1,1) = false
    options.a0 (1,1) {mustBeNumeric,mustBePositive} = 1
end

% get parameters
varfun = options.varfun;
distbins = options.distbins;
minstreams = options.minstreams;
a0 = options.a0;
    

% get node attribute list with flow accumulation values
a = ezgetnal(S,A);

% get node attribute list with elevation values
z = ezgetnal(S,DEM,'double');

% z must be double precision
z = double(z);
if options.zerobaselevel
    z = zerobaselevel(S,z);
end
copyz = z;
z = z - min(z) + 1;

if options.optimize
% calculate stream 
label = labelreach(S);

% choose optimizer
switch lower(options.funoptim)
    case 'fminsearch'
        mn = fminsearch(@(mn) chivar(mn), log(options.mn0));
        mn = exp(mn);
end
else 
    mn = options.mn0;
end

if nargout > 1 || options.plot
    c   = chitransform(S,a,'mn',mn,'a0',a0);
    [~,~,bin] = histcounts(c,distbins);
    zm  = accumarray(bin,copyz,[],@mean,0);
    zsd = accumarray(bin,copyz,[],@std,0);
    cm  = accumarray(bin,c,[],@mean,0);
    
end

if options.plot
    plotdz(S,copyz,'distance',c,'color',[.7 .7 .7])
    hold on
    errorbar(cm,zm,zsd,'k.')
    plot(cm,zm,'ks-');
    hold off
    xlabel('\chi [m]')
end
    

%% objective function
function v = chivar(mn)
    
    % calculate chitransform
    c = chitransform(S,a,'mn',exp(mn),'a0',a0);
    % bin c values into distance bins.
    [~,~,bin] = histcounts(c,distbins);
    
    % calculate the average values within each bin and tributary label
    melev = accumarray([label bin],z,[],@mean,0,true);
    
    % calculate cross-tributary variability using varfun
    n     = sum(spones(melev))';
    [~,j,melev] = find(melev);
    v    = accumarray(j(n(j)>=minstreams),melev(n(j)>=minstreams),...
                      [numel(n) 1],varfun,nan);
    
    % exclude those variances that have  less than minstreams values              
    I    = n < minstreams;
    v(I) = [];
    n(I) = [];
    
    % calculate the weighted average of all variances. The weighting is
    % according to the number of streams used to calculate the variance in
    % each bin.
    n = n/sum(n);
    v    = n'*v;
end
end
