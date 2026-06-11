function z = cumtrapz(S,G,options)

%CUMTRAPZ Cumulative trapezoidal numerical integration along a stream network
%
% Syntax
%
%     z = cumtrapz(S,G)
%     z = cumtrapz(S,g)
%     z = cumtrapz(S,g,'distance',d,'uselibtt',tf)
%
% Description
%
%     cumtrapz(S,G) computes the cumulative integral of G with respect to
%     the distance along the stream network S using trapezoidal
%     integration. S must be a stream network. The second input must either
%     be a GRIDobj G or a node attribute list g.
%     
% Input arguments
%
%     S     STREAMobj
%     G     GRIDobj
%     g     node attribute list
%
%     Parameter name/value pairs
%     
%     'distance' - {distance(S,'node_to_node')}
%          a node attribute list with inter-node distances between a node
%          and its downstream neighbor. 
%     'uselibtt' - {true},false
%          if true, then cumtrapz will be calculated using functions in
%          libtopotoolbox. If false, then native MATLAB code will be used.
%
% Output arguments
%
%     z     node attribute list
%
% Example 1
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     A = flowacc(FD);
%     a = getnal(S,A)*DEM.cellsize^2;
%     ghat = 1./(a.^0.45);
%     z = cumtrapz(S,ghat);
%     plotdz(S,z)
%
% Example 2
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     S = klargestconncomps(S);
%     A = flowacc(FD);
%     c = chitransform(S,A);
%     dc = diff(S,c);
%     ghat = smooth(S,rand(size(S.x)),'K',20);
%     ghat = normalize(ghat,'range',[0 1]);
%     ix = randlocs(S,20);
%     ghat(ismember(S.IXgrid,ix)) = 20;
%     z = cumtrapz(S,ghat,'distance',dc);
%     [~,zb] = zerobaselevel(S,DEM);
%     plotdz(S,z+zb)
%
% See also: STREAMobj, STREAMobj/gradient
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 8. June, 2026

arguments
    S   STREAMobj
    G   {mustBeGRIDobjOrNalOrScalar(G,S)}
    options.distance = distance(S,'node_to_node')
    options.uselibtt (1,1) = (true && haslibtopotoolbox)
end

% get node attribute list with elevation values
if isa(G,"GRIDobj")
    g = ezgetnal(S,G,underlyingType(G));
elseif isscalar(G)
    g = ezgetnal(S,G,class(G));
else
    g = G;
end

% get inter-node distances as edge properties
d = options.distance;
validateattributes(d,{'numeric'},{'size',size(S.x)});
d = d(S.ix);

ix  = S.ix;
ixc = S.ixc;

% Using libtt has a couple of requirements
uselibtt = options.uselibtt && ...
           haslibtopotoolbox && ...
           isfloat(g) && ...
           isreal(g);

if uselibtt
    % Use libtopotoolbox
    if isa(g,"single")
        z = tt_streamquad_trapz_f32(g,int64(ix-1),int64(ixc-1),single(d));
        return
    elseif isa(g,"double")
        z = tt_streamquad_trapz_f64(g,int64(ix-1),int64(ixc-1),single(d));
        return
    end
else
    % Use MATLAB code
    z = zeros(size(g));
    g = double(g);
    for r = numel(ix):-1:1
        z(ix(r)) = z(ixc(r)) + (g(ixc(r))+(g(ix(r))-g(ixc(r)))/2)*d(r);
    end

end

%% Is cumtrapz really working correctly. Here is the test:
% d  = S.distance;
% d2 = cumtrapz(S,ones(size(S.x)));
% tf = isequal(d,d2)


