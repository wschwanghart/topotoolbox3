function varargout = getoutlets(FD)

%GETOUTLETS Retrieve outlets of a flow network (FLOWobj)
%
% Syntax
%
%     ix = getoutlets(FD)
%     [x,y] = getoutlets(FD)
%
% Description
%
%     This function takes a FLOWobj as input and identifies the outlet
%     locations of the drainage network. Outlets are defined as nodes where
%     flow exits the domain or where no downstream neighbor exists within
%     the flow network. The function analyzes the flow direction and
%     connectivity stored in the FLOWobj to detect these terminal points
%     and returns their positions either as linear indices (ix) referencing
%     the DEM grid cells or as coordinate vectors (x, y) corresponding to
%     their spatial locations. This enables users to extract, visualize, or
%     further analyze watershed outlets within a hydrological modeling
%     workflow.
%
% Input arguments
%
%     FD   FLOWobj
%
% Output arguments
%
%     ix   linear index into DEM from which FD was derived
%     x,y  coordinate vectors
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     [x,y] = getoutlets(FD);
%     A  = flowacc(FD);
%     imageschs(DEM,sqrt(A),'colormap',flowcolor);
%     hold on
%     plot(x,y,'ok','MarkerFaceColor','k')
%     hold off
%     padextent(500);
%
% See also: FLOWobj/drainagebasins, STREAMobj/streampoi
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 8. October, 2025

arguments
    FD FLOWobj
end

nargoutchk(1,2)

ix  = FD.ix;
ixc = FD.ixc;
I   = false(FD.size);

% Outlets are part of the drainage network but are not members of the
% giving nodes.
I(ixc) = true;
I(ix) = false;

if nargout == 1
    varargout{1} = find(I);
elseif nargout == 2
    [varargout{1},varargout{2}] = findcoord(GRIDobj(FD,I));
end
