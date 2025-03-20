function [a,D] = hillslopearea(S,FD,W)

%HILLSLOPEAREA Upslope hillslope area for each stream pixel 
%
% Syntax
%
%     a = hillslopearea(S,FD)
%     [a,D] = hillslopearea(S,FD)
%     ... = hillsloparea(S,FD,W)
%     
% Description
%
%     hillslopearea returns the upslope hillslope area of each river pixel
%     in S. Compared to flow accumulation, this function stops accumulation
%     along river pixels so that the accumulated flow calculated for each
%     river pixels includes only the hillslope pixels but not those further
%     upstream along the stream network. The function also returns the
%     GRIDobj D which contains labels that for each channel site and the
%     hillslope area that drains into this site without passing another
%     upstream channel site (see also Fig. 4 in Hergarten 2021).
%
% Input arguments
%
%     S     STREAMobj
%     FD    FLOWobj
%     W     GRIDobj used to weigh upstream areas (see flowacc)
%
% Output arguments
%
%     a     node-attribute list with hillslope areas
%     D     GRIDobj with drainage basins for each river pixel
%
% Reference: 
% 
%     Hergarten, S.: Rivers as linear elements in landform evolution
%     models, Earth Surface Dynamics, 8, 367–377,
%     https://doi.org/10.5194/esurf-8-367-2020, 2020.
%
% See also: FLOWobj/flowacc, FLOWobj/upslopestats
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 21. February, 2025    

arguments
    S    STREAMobj
    FD   FLOWobj
    W    GRIDobj = GRIDobj(FD)+1;
end

A = upslopestats(FD,W,'sum',S);
a = getnal(S,A);

if nargout == 2
    D = drainagebasins(FD,S.IXgrid);
end