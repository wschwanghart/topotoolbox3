function z = lowerenv(S,z,ix,options)

%LOWERENV Lower envelope of a channel length profile
%
% Syntax
%
%     zl = lowerenv(S,z)
%     zl = lowerenv(S,z,kn)
%     zl = lowerenv(___, uselibtt = false)
%
% Description
%
%     lowerenv returns the lower envelope, i.e. the lower convex hull of a
%     length profile given by the stream network S and elevation z.
%
% Input arguments
%
%     S      STREAMobj
%     z      elevation (node attribute list) or GRIDobj. Ideally, the
%            elevations should be preprocessed using imposemin, quantcarve
%            or crs. Otherwise, the returned elevations may not be strictly
%            downward decreasing.
%     kn     logical vector (node attribute list) with knickpoints. Default
%            is [].
%
%     Parameter name/value pairs
%
%     uselibtt    true or {false}. If true, the function will use
%                 libtopotoolbox (if available).     
%
% Output arguments
%
%     zl     node-attribute list with elevations (double)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD    = FLOWobj(DEM);
%     S = klargestconncomps(STREAMobj(FD,'minarea',1000));
%     S = trunk(S);
%     z = imposemin(S,DEM);
%     zs = lowerenv(S,z);
%     plotdz(S,zs);
%     hold on
%     plotdz(S,DEM)
%
%     % Example with knickpoints
%     [~,kp] = knickpointfinder(S,DEM,'split',false,...
%                                     'tol',20,'plot',false,...
%                                     'verbose',false);
%     zs = lowerenv(S,z,kp.nal);
%     plotdz(S,zs);
%     hold on
%     plotdz(S,DEM)    
%     
% See also: STREAMobj/knickpointfinder, STREAMobj/imposemin, 
%           STREAMobj/quantcarve, STREAMobj/crs    
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 23. July, 2025

arguments
    S STREAMobj
    z {mustBeGRIDobjOrNal(z,S)}
    ix = []
    options.uselibtt (1,1) = false
end

z = ezgetnal(S,z);

kn = false(size(S.x));
if nargin == 3
    if isempty(ix)
        % do nothing
    elseif islogical(ix)
        kn = ix;
    else
        kn(ix) = true;
    end
end

% nals
d = distance(S);

if options.uselibtt && haslibtopotoolbox
    z = single(z);
    z = tt_lowerenv(z, uint8(kn), single(d), ...
        int64(S.ix - 1), int64(S.ixc - 1));
    return;
end

trib = streampoi(S,'confl','logical');

nrc = numel(S.x);
ix  = S.ix;
ixc = S.ixc;

ixcix  = zeros(nrc,1);
ixcix(ix) = 1:numel(ix);

onenvelope = true(nrc,1);

for r = numel(S.ixc):-1:1
    s  = ix(r);
    ss = ixc(r);
    
    if onenvelope(s) || trib(ss)
    IX = allpred(ix,ixc,s,nrc,kn,r);
    s  = ss;
    if isempty(IX)
        continue
    end
    gg = z(IX)-z(s);
    dd = d(IX)-d(s);
    gg = gg./dd;
    [~,ii] = sortrows([gg dd],[1 -2]); 
    IX = IX(ii(1));
    g  = gg(ii(1));
    
    t = IX;
    ixcix(ss) = 0;
    while ixcix(t) ~= 0
        t2 = ixc(ixcix(t));
        z(t2) = z(t)-g.*(d(t)-d(t2));
        onenvelope(t2) = false;
        ixcix(t) = 0;
        t  = t2;
        
    end
    end
    
end

end


function IX = allpred(ix,ixc,s,nr,kn,startix)
I = false(nr,1);
I(s) = true;
for r = startix:-1:1
    I(ix(r)) = I(ix(r)) | (I(ixc(r)) & ~kn(ixc(r))) ;
end
IX = find(I);
end

