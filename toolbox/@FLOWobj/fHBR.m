function H = fHBR(FD,DEM,S,options)

%FHBR Fuzzy Height Below Ridge (fHBR)
%
% Syntax
%
%     H = fHBR(FD,DEM)
%     H = fHBR(FD,DEM,S)
%
% Description
%
%     What is the height difference between a pixel and the nearest ridge?
%     The difficulty in answering this simple question is that a ridge
%     elevation cannot be unambiguously be assocatiated to a pixel. The
%     function fHBR tries to account for this fuzziness by using multiple
%     flow directions stored in the FLOWobj FD. By using multiple flow
%     directions, a pixel potentially has multiple upstream ridge-pixels
%     that are identified as those that have no upstream neighboring
%     pixels. fHBR calculates the weighted average of the elevation of 
%     these ridge pixels.
%
% Input arguments
%
%     FD    FLOWobj (multiple or dinf)
%     DEM   digital elevation model (GRIDobj)
%     S     STREAMobj. If provided, computations will stop at stream pixels
%
%     Parameter name/value pairs
%
%     'seed'        GRIDobj or linear index for seed locations from which 
%                   distances will be calculated 
%
% Output arguments
%
%     H     Height below ridge (GRIDobj)
%
% Example: 
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'dinf');
%     [FD,S] = multi2single(FD,"minarea",500);
%     D = meanflowdist(FD,S,"direction","downstream");
%     H = fHBR(FD,DEM,S);
%     plot(D.Z(:),-H.Z(:),'.')
% 
%     % Also pick a more homogeneous area with similar slopes using
%     % createmask and plot the relation between D and H again
%     BW = createmask(DEM);
%     plot(D.Z(BW.Z),-H.Z(BW.Z),'.')
%
% See also: FLOWobj, FLOWobj/meanflowdist, FLOWobj/fHAND,
%           FLOWobj/slopeposition
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 2. September, 2026

arguments
    FD   FLOWobj
    DEM  GRIDobj
    S = []
    options.seed = []
end

if ~isempty(S)
    S = STREAMobj2GRIDobj(S);
    I = S.Z(FD.ix) & S.Z(FD.ixc);
    I = ~I;
    FD.ix = FD.ix(I);
    FD.ixc = FD.ixc(I);
    if ismulti(FD)
        FD.fraction = FD.fraction(I);
    end
end

if ~isempty(options.seed)
    I = influencemap(FD,options.seed);
    I = I.Z(FD.ix);
    FD.ix = FD.ix(I);
    FD.ixc = FD.ixc(I);
    if ~isempty(FD.fraction)
        FD.fraction = FD.fraction(I);
    end
    clear I
end

ix = FD.ix;
ixc = FD.ixc;
fraction = FD.fraction;

if isempty(fraction)
    fraction = ones(size(ix));
end

% Identify pixels with no upstream neighbors
I = false(DEM.size);
I(ix) = true;
I(ixc) = false;

Z = DEM.Z;

Z(~I) = 0;
A     = flowacc(FD).Z;
w     = zeros(size(ix));
wsum  = zeros(DEM.size);

for r = 1:numel(FD.ix)
    w(r) = A(ix(r)) * fraction(r);
    wsum(ixc(r)) = wsum(ixc(r)) + w(r);
end

for r = 1:numel(ix)
    Z(ixc(r)) = Z(ixc(r))+w(r)*Z(ix(r))/wsum(ixc(r));
end

H = GRIDobj(DEM,Z)-DEM;
H = max(H,0);

