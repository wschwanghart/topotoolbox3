function FD = clip(FD,MASK)

%CLIP Clip FLOWobj
%
% Syntax
%
%     FDc = clip(FD,MASK)
%
% Description
%
%     Some analysis will likely not require the entire flow network. Hence,
%     the function clip enables you to clip a flow network stored in FD by
%     removing all links with start and endpoints outside the mask MASK
%     provided as GRIDobj.
%
% Input arguments
%
%     FD       FLOWobj
%     MASK     GRIDobj (logical)
%     
% Output arguments
%
%     FDc   clipped FLOWobj
%
% See also: FLOWobj2cell
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 10. March, 2025 

arguments
    FD  FLOWobj
    MASK   GRIDobj {validatealignment(FD,MASK)} 
end

validatealignment(FD,MASK)
II = MASK.Z(FD.ix) & MASK.Z(FD.ixc);
FD.ix = FD.ix(II);
FD.ixc = FD.ixc(II);
if ismulti(FD)
    FD.fraction = FD.fraction(II);
end
