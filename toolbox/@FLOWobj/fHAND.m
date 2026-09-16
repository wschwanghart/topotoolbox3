function V = fHAND(FD,DEM,S,options)

%fHAND Fuzzy height above nearest drainage (HAND)
%
% Syntax
%
%     V = fHAND(FD,DEM,S)
%     V = fHAND(FD,DEM,S,pn,pv,...)
%
% Description
%
%     fHAND is a fuzzy implementation of the height above nearest drainage 
%     algorithm (HAND). Relying on multiple flow directions stored in the
%     FLOWobj FD, every hillslope pixel is associated with one or more
%     drainage pixels stored in the STREAMobj S. The resulting GRIDobj V
%     contains the weighted average of the difference between the hillslope
%     pixel's elevation and that of the downslope drainage pixels. 
%
% Input arguments
%
%     FD    FLOWobj (multi or single)
%     S     STREAMobj
%     DEM   digital elevation model (GRIDobj)
%
%     Parameter name/value pairs
%
%     'calcdiff'  {true} or false. If false, the algorithm returns the
%                 average drainage pixel elevation for each hillslope 
%                 pixel.
%     'setstreamsnan' {false} or true. If true, stream pixels in V are set
%                 to nan.
%
% Output arguments
%
%     V     GRIDobj with (average) heights above nearest drainage 
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'multi');
%     [FD,S] = multi2single(FD,'minarea',1000);
%     V = fHAND(FD,DEM,S);
%     imageschs(DEM,V)
%
% Reference: Roering JJ et al. 2026. Shallow Landslides Align With
%          Atmospheric Rivers in Coastal Steeplands. Geophysical Research 
%          Letters 53, e2026GL124294. DOI: 10.1029/2026GL124294
%
% See also: FLOWobj, FLOWobj/meanflowdist, FLOWobj/fHBR,
%           FLOWobj/slopeposition
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 2. September, 2026

arguments
    FD   FLOWobj
    DEM  GRIDobj
    S    = []
    options.calcdiff = true
    options.setstreamsnan = false
end

if isa(S,"STREAMobj")

    DEMS = STREAMobj2GRIDobj(S,DEM);
    DEMS.Z(isnan(DEMS.Z)) = 0;
    ISSTREAM = STREAMobj2GRIDobj(S);
    streamlocs = S.IXgrid;
else
   % get outlet locations
   outlets = true(DEM.size);
   for r = numel(FD.ix):-1:1
       outlets(FD.ix(r)) = false;
   end
   streamlocs = find(outlets);
   DEMS = GRIDobj(DEM,underlyingType(DEM));
   DEMS.Z(streamlocs) = DEM.Z(streamlocs);
   ISSTREAM = GRIDobj(DEM,outlets);
end


% There are some pixels that do not have downstream stream pixels. Usually,
% these are part of smaller catchment along the DEM edges. We remove links
% from the flow network that connect any of these pixels.
I   = dependencemap(FD,streamlocs);
ii   = I.Z(FD.ix) & I.Z(FD.ixc);
FD.ix = FD.ix(ii);
FD.ixc = FD.ixc(ii);
if isempty(FD.fraction)
    FD.fraction = ones(size(FD.ix));
else
    FD.fraction = FD.fraction(ii);
end

% Accumulated fractions 
FR  = zeros(DEMS.size);
A   = flowacc(FD).Z;
A   = A(FD.ix);

for r = numel(FD.ix):-1:1
    % FR(FD.ix(r)) = FD.fraction(r) + FR(FD.ix(r));
    FR(FD.ix(r)) = A(r).*FD.fraction(r) + FR(FD.ix(r));
end

for r = numel(FD.ix):-1:1
    if ISSTREAM.Z(FD.ix(r))
        continue
    end

    % DEMS.Z(FD.ix(r)) = FD.fraction(r).*DEMS.Z(FD.ixc(r)) .* 1/FR(FD.ix(r)) ...
    %                     + DEMS.Z(FD.ix(r));
        
    DEMS.Z(FD.ix(r)) = A(r).*FD.fraction(r).*DEMS.Z(FD.ixc(r)) .* 1/FR(FD.ix(r)) ...
                        + DEMS.Z(FD.ix(r));
end

if options.calcdiff
    V = DEM-DEMS;
else
    V = DEMS;
end

if options.setstreamsnan
V = clip(V,~ISSTREAM);
end
V.Z(~I.Z) = nan;