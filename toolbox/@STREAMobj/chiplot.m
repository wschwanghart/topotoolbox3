function OUT = chiplot(S,DEM,A,options)

%CHIPLOT Chi analysis for bedrock river analysis
%
% Syntax
%
%     C = chiplot(S,DEM,A)
%     C = chiplot(S,z,a)
%     C = chiplot(S,DEM,A,pn,pv,...)
%
% Description
%
%     CHI plots are an alternative to slope-area plots for bedrock river
%     analysis. CHI plots are based on a transformation of the horizontal
%     coordinate that converts a steady-state river profile into a straight
%     line with a slope that is related to the ratio of the uplift rate to
%     the erodibility (Perron and Royden 2012).
%
% Input arguments
%
%     S     instance of STREAMobj. The stream network must consist of only 
%           one connected component (only one outlet may exist!)
%     DEM   digital elevation model (GRIDobj)
%     A     flow accumulation as calculated by flowacc (GRIDobj)
%     z     elevation values for each node in S (node attribute list)
%     a     flowaccumulation values for each node in S (node attribute list)
%
%     Parameter name/value pairs {default}
%
%     'a0': {1e6}
%     reference area in m^2
%
%     'mn':  {[]}, scalar      
%     mn is the ratio of m and n in the stream power equation. The value
%     ranges usually for bedrock rivers between 0.1 and 0.5. If empty, it
%     is automatically found by a least squares approach.
%
%     'mnoptim': {'fminsearch'}, 'nlinfit'
%     chiplot finds an optimal value of mn that linearizes the relation
%     between chi and elevation. By default, fminsearch is used. nlinfit
%     requires the statistics and machine learning toolbox but additionally
%     returns 95% confidence intervals of mn. These bounds must, however,
%     be taken with care since their determination violates the assumption
%     of an independent and identically distributed (i.i.d.) variable. Note 
%     also that the standard deviation of beta (betase) is determined
%     independently from mn.
%     
%     'color': color of river profiles, e.g. 'b' or [.4 .4 .6]
%
%     'trunkstream': {[]}, STREAMobj
%     instance of STREAMobj that must be a subset of S, e.g. the main river
%     in the network S. The main trunk is highlighted in the plot and can
%     be used to fit the mn ratio (see pn/pv pair 'fitto'). Note that the
%     trunkstream must end at the river network's root (outlet).
%
%     'colorMainTrunk': color of the trunk stream, e.g. 'b' or [.4 .4 .6].
%     Applies only if 'trunkstream' (see above) is supplied.
%
%     'fitto': {'all'},'ts', 
%     choose which data should be used for fitting the mn ratio.
%     'all' fits mn to all streams in S
%     'ts' fits mn only to the trunkstream which must be provided with the
%     pn-pv pair trunkstream.
%
%     'plot': {true}, false
%     plot the CHIplot.
%
%     'mnplot': {false}, true
%     plot data for various values of mn (see mnvalues)
%
%     'mnvalues': [0.2:0.1:0.8] 
%     different values of mnvalues for chiplots if mnplot is set to true.
%
%     'normchi': {false} or true
%     if true, then the horizontal distance in the mnplot will be
%     normalized to range between 0 and 1.
%
% Output arguments
%
%     C     structure array that contains
%     .mn       ratio of m and n
%     .mnci     95% intervals of mn (empty if 'mnoptim' is set to
%               fminsearch)
%     .beta     slope of the best fit line
%     .betase   standard error of beta
%     .a0       reference area
%     .ks       channel steepness index
%     .chi      CHI values 
%     .elev     elevation
%     .elevbl   elevation above baselevel
%     .distance distance from outlet
%     .pred     predicted elevation
%     .res      residual elevation
%
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','carve');
%     S   = STREAMobj(FD,'minarea',1000);
%     S   = klargestconncomps(S);
%     c   = chiplot(S,DEM,flowacc(FD));
%
% See also: flowpathapp, STREAMobj, FLOWobj/flowacc, STREAMobj/trunk,
%           STREAMobj/modify, STREAMobj/getnal, STREAMobj/chitransform
%
% References:
%     
%     Perron, J. & Royden, L. (2013): An integral approach to bedrock river 
%     profile analysis. Earth Surface Processes and Landforms, 38, 570-576.
%     [DOI: 10.1002/esp.3302]
%     
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de) and 
% %       Karina Marques (https://github.com/karina-marques)
% Date: 12. November, 2024

% update 11. June, 2014
% supports node attribute lists as input data (DEM,A)
% update 4. October, 2016
% lets you choose the algorithm to find the optimal value of mn. 
% update 30. June, 2022
% added options to color the output by Karina Marques 

arguments
    S  STREAMobj
    DEM {mustBeGRIDobjOrNal(DEM,S)}
    A {mustBeGRIDobjOrNal(A,S)}
    options.color = []
    options.colorMainTrunk = []
    options.mn = []
    options.trunkstream = []
    options.plot (1,1) = true
    options.mnplot (1,1) = false
    options.mnvalues (1,:) {mustBePositive,mustBeNonempty} = [0.2:0.1:0.8]
    options.fitto {mustBeMember(options.fitto,{'all','ts'})} = 'all'
    options.mnoptim {mustBeMember(options.mnoptim,...
        {'fminsearch','nlinfit'})}= 'fminsearch'
    options.a0 (1,1) {mustBePositive} = 1e6
    options.betamethod {mustBeMember(options.betamethod,{'ls','lad'})} = 'ls'
    options.mnmethod {mustBeMember(options.mnmethod,{'ls','lad'})} = 'ls'
    options.normchi (1,1) = false
end

fitto = options.fitto;
color = options.color;
colorMainTrunk = options.colorMainTrunk;
betamethod = validatestring(options.betamethod,{'ls','lad'});
mnmethod   = validatestring(options.mnmethod,{'ls','lad'});

if options.plot
    ax = gca;
end

% Make sure that a trunk stream is supplied if 'ts' is chosen
switch fitto
    case 'ts'
        if isempty(options.trunkstream)
            error('TopoToolbox:wronginput',...
                ['You must supply a trunkstream, if you use the parameter \n'...
                'fitto together with the option ts']);
        end
end

% nr of nodes in the entire stream network
outlet = streampoi(S,'outlet','logical');
if nnz(outlet)>1
    % there must not be more than one outlet (constraint could be removed
    % in the future).
    error('TopoToolbox:chiplot',...
        'The stream network must not have more than one outlet');
end

% reference drainage area
a0   = options.a0; % m^2
% elevation values at nodes
zx   = ezgetnal(S,DEM,'double');
zb   = zx(outlet);

a = ezgetnal(S,A)*(A.cellsize.^2); 
a = a0./a;

% x is the cumulative horizontal distance in upstream direction
x    = S.distance;

% use trunkstream for fitting and display
switch fitto
    case 'all'
        SFIT = S;
        Lib  = true(size(x));

    case 'ts'
        SFIT = options.trunkstream;
        [Lia,Lib] = ismember(SFIT.IXgrid,S.IXgrid);
        if any(~Lia)
            error('TopoToolbox:chiplot',...
                ['The main trunk stream must be a subset of the stream network.\n'...
                 'Map a trunk stream with flowpathtool and use the STREAMobj as \n'...
                 '3rd input argument ( flowpathtool(FD,DEM,S) ).']) ;
        end
end

% find values of the ratio of m and n that generate a linear Chi plot
% uses fminsearch
if isempty( options.mn )
    mn0  = 0.5; % initial value
    % fminsearch is a nonlinear optimization procedure that doesn't require
    % the statistics toolbox
    switch options.mnoptim
        case 'fminsearch'
            mn   = fminsearch(@mnfit,mn0);
            ci   = [];
        case 'nlinfit'
            ztest   = zx(Lib)-zb;
            ztest   = ztest./max(ztest);
            [mn0,R,J] = nlinfit(a,ztest,@mnfit,mn0);
            ci = nlparci(mn0,R,'jacobian',J);
    end
else
    % or use predefined mn ratio.
    mn   = options.mn;
    ci   = [];
end

% plot different values of mn
if options.mnplot
    mntest = options.mnvalues;
    cvec = jet(numel(mntest));
    figure('DefaultAxesColorOrder',cvec);
    
    for r = 1:numel(mntest)
        c = chitransform(S,a/S.cellsize,"a0",options.a0,"mn",mntest(r));
        if options.normchi
            c = c/max(c);
        end
        plotdz(S,zx,'distance',c,'color',cvec(r,:))
        hold on
    end
    hold off
    
    if options.normchi
        xlabel('\chi [m] (normalized)')
    else
        xlabel('\chi [m]')
    end
    ylabel('Elevation [m]');
    title('\chi plots for different values of mn')
    legnames = cellfun(@(x) num2str(x),num2cell(mntest),'uniformoutput',false);
    legend(legnames);
end

% calculate chi
chi = cumtrapz(S,a.^mn);
% now use chi to fit beta
switch betamethod
    case 'ls'
        % least squares
        beta = chi(Lib)\(zx(Lib)-zb);
    case 'lad'
        % least absolute deviations
        beta = fminsearch(@(b) sum(abs(b*chi(Lib) - (zx(Lib)-zb))),0.0334);
end

n    = nnz(Lib);
SSE  = sum((chi(Lib)*beta - (zx(Lib)-zb)).^2);
SSZ  = sum((zx(Lib)-mean(zx(Lib))).^2);
R2   = 1-(SSE/SSZ);

betase = sqrt((SSE/(n-2))./(sum((chi(Lib)-mean(chi(Lib))).^2)));


if options.plot
    % plot results
    order = S.orderednanlist;
    I     = ~isnan(order);
    c     = nan(size(order));
    c(I)  = chi(order(I));
    zz    = nan(size(order));
    zz(I) = zx(order(I));

    if isempty(color)

        colororderindex = mod(ax.ColorOrderIndex, size(ax.ColorOrder,1));
        if colororderindex==0; colororderindex=size(ax.ColorOrder,1); end
        color = ax.ColorOrder(colororderindex,:);
    end
    
    plot(ax,c,zz,'-','color', color);
    hold on

    if ~isempty( options.trunkstream )
        switch options.fitto
            case 'all'
                % check trunkstream
        
                ST    = options.trunkstream;
                [Lia,Lib] = ismember(ST.IXgrid,S.IXgrid);
        
                if any(~Lia)
                    error('TopoToolbox:chiplot',...
                        ['The main trunk stream must be a subset of the stream network.\n'...
                        'Map a trunk stream with flowpathtool and use the STREAMobj as \n'...
                        '3rd input argument ( flowpathtool(FD,DEM,S) ).']) ;
                end
            case 'ts'
                 ST = SFIT;
        end
        
        order = ST.orderednanlist; 
        I     = ~isnan(order);
        c     = nan(size(order));
        chifit = chi(Lib);
        zxfit  = zx(Lib);
        c(I)  = chifit(order(I));
        zz    = nan(size(order));
        zz(I) = zxfit(order(I));

        if isempty(colorMainTrunk)
            colorMainTrunk = color;
        end

        plot(ax,c,zz,'color', colorMainTrunk,'LineWidth',2);
    end
    
    refline(ax,beta,zb);
    hold off
    xlabel('\chi [m]')
    ylabel('elevation [m]');
end

% write to output array
if nargout == 1
    
    OUT.mn   = mn;
    OUT.mnci = ci;
    OUT.beta = beta;
    OUT.betase = betase;
    OUT.a0   = a0;
    OUT.ks   = beta*a0^mn;
    OUT.R2   = R2;
    OUT.x_nal    = S.x;
    OUT.y_nal    = S.y;
    OUT.chi_nal  = chi;
    OUT.d_nal    = S.distance;
    
    [OUT.x,...
     OUT.y,...
     OUT.chi,...
     OUT.elev,...
     OUT.elevbl,...
     OUT.distance,...
     OUT.pred,...
     OUT.area] = STREAMobj2XY(S,chi,DEM,zx-zb,S.distance,beta*chi,A.*(S.cellsize^2));
     OUT.res   = OUT.elevbl-OUT.pred;

end
    


%% fitting function
function sqres = mnfit(varargin)
    
    if nargin == 1
        mn = varargin{1};
    elseif nargin == 2
        mn = varargin{1};
        a  = varargin{2};
    end

% calculate chi with a given mn ratio
% and integrate in upstream direction
CHI = cumtrapz(SFIT,a(Lib).^mn);
% normalize both variables
CHI = CHI ./ max(CHI);

if nargin == 2
    sqres = CHI;
    return
end

z   = zx(Lib)-zb;
z   = z./max(z);
% calculate the residuals and minimize their squared sums

switch mnmethod
    case 'ls'
        sqres = sum((CHI - z).^2);
    case 'lad'
        sqres = sum(sqrt(abs(CHI-z)));
end 

end

end






