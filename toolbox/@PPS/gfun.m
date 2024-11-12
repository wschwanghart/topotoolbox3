function g = gfun(P,options)

%GFUN G-function (nearest inter-point distance distribution)
%
% Syntax
%
%     g = gfun(P)
%     d = gfun(P,pn,pv,...)
%
% Description
%
%     gfun calculates the cumulative probability distribution of nearest
%     neighbor distances and compares it to the distribution of a
%     homogeneous or inhomogeneous Poisson process.
%
% Input arguments
%
%     P       instance of PPS
%
%     Parameter name/value pairs
%
%     'nsim'    number of bootstrap simulations
%     'model'   mdl (see PPS/fitloglinear) or density (see PPS/density)
%     'plot'    {false} or true. 'plot' is true, if the function is called
%               without output arguments
%
% Output arguments
%
%     d       matrix or vector of distances
%
% Example
%
%
%
%
% See also: PPS
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 12. November, 2024

arguments
    P PPS
    options.nsim (1,1) {mustBeNonnegative} = 0
    options.model = []
    options.plot (1,1) = false
end

plotit = nargout == 0 || options.plot;

d      = pointdistances(P,'output','matrix');
d(eye(size(d))>0) = inf;
d      = min(d,[],2);
d      = sort(d,'ascend');

nsim   = options.nsim;
if nsim > 0
    PPSIM = cell(nsim,1);
    for r = 1:nsim
        if isempty(options.model)
            PSIM = random(P);
        else
            if isa(options.model,'GeneralizedLinearModel')
                PSIM = random(P,options.model);
            else
                PSIM = simulate(P,'intensity',options.model);
            end
        end
        PPSIM{r} = PSIM.PP;
            
    end
    G  = as(P,'graph');
    D  = cell(nsim,1);
    maxD = cell(nsim,1);
    
    % Parallel waitbar
    DQ = parallel.pool.DataQueue;
    h  = waitbar(0, 'Please wait ...');
    afterEach(DQ, @nUpdateWaitbar);
    
    counter = 1;
    
    parfor r = 1:nsim
        % point pattern
        pp = PPSIM{r};
        % distances on graphs do not allow multiple points at one location
        [pp,~,locb] = unique(pp,'stable');
        % calculate distances
        dd = distances(G,pp,pp,'Method','positive');
        % undo unique
        dd  = dd(locb,locb);
        % set distances on the main diagonal to inf
        dd(eye(size(dd))>0) = inf;
        % calculate minimum distances
        dd  = min(dd,[],2);
        % sort distances 
        dd  = sort(dd,'ascend');
        % ... and place them in the cell array 
        D{r} = dd;
        maxD{r} = dd(end);
        
        % Send to waitbar
        send(DQ, r);
    end
    
    % Close waitbar
    close(h)
    
    mind = 0;
    maxd = vertcat(maxD{:});
    maxd = max(maxd,max(d));
    
    if plotit
    for r = 1:numel(D)
        plot(D{r},(1:numel(D{r}))/numel(D{r}),'color',[.8 .8 .8]);
        hold on
    end
    end
    
    
end
if plotit
plot(d,(1:numel(d))/numel(d),'k');
end

g = [];

    function nUpdateWaitbar(~)
        waitbar(counter/nsim, h);
        counter = counter + 1;
    end
end