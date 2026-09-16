function SP = slopeposition(FD,S,DEM)

%SLOPEPOSITION Slope position index (relative relief)
%
% Syntax
%
%     SP = slopeposition(FD,S)
%     SP = slopeposition(FD,S,DEM)
%
% Description
%
%     This function computes the slope position index for each pixel in a
%     terrain model, returning a normalized value between 0 and 1 that
%     represents the relative position along a hillslope. Values close to 0
%     indicate pixels located on toe slopes and valley bottoms, while
%     values close to 1 represent upper slopes or ridges. Intermediate
%     values describe mid-slope positions, enabling a continuous
%     representation of terrain position rather than discrete classes.
%
%     The resulting slope position map provides a standardized and
%     spatially consistent way to characterize topographic structure and
%     landscape organization. It can be used in geomorphological analysis,
%     landslide susceptibility analysis, hydrological modeling, soil and
%     erosion studies, or ecological applications where distinguishing
%     between upper, middle, and lower slope zones is important.
%
%     slopeposition(FD,S) derives slope position based on the horizontal
%     distance along the slope by comparing mean upstream path length above
%     each pixel (Du) with the mean downstream path length to the nearest
%     drainage below each pixel (Dd). SP is then calculated by Du/(Du+Dd).
%
%     slopeposition(FD,S,DEM) derives slope position based on the digital
%     elevation model DEM. The function calculates for each pixel the
%     height below ridge (Hm), a metric based on multiple flow directions.
%     In addition, it calculates fHAND, a fuzzy version of the HAND (height
%     above nearest drainage) based on multiple flow directions (Hp). SP
%     is then calculated by Hp/(Hp+Hm).
%
% Input arguments
%
%     FD     Multiple flow direction object (FLOWobj), preferably generated
%            using multi2single with an area threshold (see example)
%     S      STREAMobj as returned by multi2single
%     DEM    optional DEM
%
% Output arguments
%
%     SP     GRIDobj of slope position index
%
% Example 1: Slope position based on distance
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'multi');
%     [FD,S] = multi2single(FD,'minarea',200);
%     SP  = slopeposition(FD,S);
%     imageschs(DEM,SP)
%     hold on
%     plot(S,'w')
%     [x,y] = contour(SP,5);
%     plot(x,y,'k')
%
% Example 2: Slope position based on elevation
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'multi');
%     [FD,S] = multi2single(FD,'minarea',200);
%     DEM = imposemin(S,DEM);
%     SP  = slopeposition(FD,S,DEM);
%     imageschs(DEM,SP)
%     hold on
%     plot(S,'w')
%     [x,y] = contour(SP,5);
%     plot(x,y,'k')
%
% Reference: Roering JJ et al. 2026. Shallow Landslides Align With
%          Atmospheric Rivers in Coastal Steeplands. Geophysical Research 
%          Letters 53, e2026GL124294. DOI: 10.1029/2026GL124294
%
% See also: FLOWobj/meanflowdist, FLOWobj/fHAND, FLOWobj/fHBR
%           
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 2. September, 2026

arguments
    FD   FLOWobj
    S    STREAMobj
    DEM  {mustBeGRIDobjOrEmpty} = []
end


if isempty(DEM)
    % Mean flow distance in downstream direction
    Dd = meanflowdist(FD,S);

    % Mean flow distance in upstream direction
    Du = meanflowdist(FD,S,"direction","upstream","A",flowacc(FD));

    % Slope position is here defined as Dd/(Dd+Du)
    SP = Du/(Dd+Du);
else
    % Calculate height below ridge (hbr)
    Hp = fHBR(FD,DEM,S);

    % Calculate fuzzy height above nearest drainage (HAND)
    Hm = fHAND(FD,DEM,S);

    % Slope position is here defined as SP = H/(H+H2)
    SP = Hm/(Hp+Hm);
end

SP.Z(S.IXgrid) = 0;

