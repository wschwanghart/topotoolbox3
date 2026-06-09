function z = cumtrapz(S,G,options)

%CUMTRAPZ Cumulative trapezoidal numerical integration along a stream network
%
% Syntax
%
%     z = cumtrapz(S,G)
%     z = cumtrapz(S,g)
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
%     'uselibtt' - {true},false
%          if true, then cumtrapz will be calculated using functions in
%          libtopotoolbox. If false, then native MATLAB code will be used.
%
% Output arguments
%
%     z     node attribute list
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     S = STREAMobj(FD,'minarea',1000);
%     A = flowacc(FD);
%     a = getnal(S,A)*DEM.cellsize^2;
%     ghat = 1./(a.^0.45);
%     z = cumtrapz(S,ghat);
%     plotdz(S,z)
%
% See also: STREAMobj, STREAMobj/gradient
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 8. June, 2026

arguments
    S   STREAMobj
    G   {mustBeGRIDobjOrNal(G,S)}
    options.uselibtt (1,1) = (true && haslibtopotoolbox)
end

% get node attribute list with elevation values
if isa(G,"GRIDobj")
    g = ezgetnal(S,G,underlyingType(G));
else
    g = G;
end

% get inter-node distances as edge properties
d = distance(S,'node_to_node');
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


