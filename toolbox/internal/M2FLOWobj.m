function FD = M2FLOWobj(M,wf,options)

%M2FLOWobj Convert transfer matrix to FLOWobj
%
% Syntax
%
%     FD = M2FLOWobj(M,wf)
%     FD = M2FLOWobj(M,georef)
%
% Description
%
%     The function converts a sparse transfer matrix M to a FLOWobj. The 
%     transfer matrix was used in TopoToolbox 1 to represent flow
%     directions. M is sparse and has as many rows and columns as there are
%     pixels in the DEM from which it was derived. Nonzero elements in M
%     represent fractions that are transferred from a giver (or donor)
%     pixel to receiver pixels. Giver pixels are along rows, and receiver
%     pixels along columns. If M(5,7) = 1, for example, the pixel with the
%     linear index 5 in the DEM, transfer all its water/sediment to its
%     downstream neighbor with the linear index 7.
%     
%     Typically, values along rows add up to one. Otherwise, a donor pixel
%     transfers more or less to its downstream neighbors than it receives
%     and/or contains.
%
% Input arguments
%
%     M         sparse transfer matrix
%     wf        world-file matrix (see GRIDobj.wf)
%     georef    map or geo referencing object (e.g. MapCellsReference)
%     
% Output arguments
%
%     FD        FLOWobj
%
% See also: FLOWobj
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 6. September, 2024

arguments
    M    {mustBeSparse}
    wf   
    options.size = estimateSizeFromM(M)
end

FD = FLOWobj;
if isequal(size(wf),[2 3])
    FD.wf = wf;
    FD.georef = [];
else
    FD.wf = worldFileMatrix(wf);
    FD.georef = wf;
end

FD.cellsize = wf2cellsize(FD.wf);
FD.size     = options.size;

% topologically sort M and generate list of ordered
% vertex links

% Is it a multiple flow direction algorithm?
mtf = any(sum(spones(M),2)>1);

if mtf
    [ix,ixc,fr] = find(M);

    FD.fraction = fr;
    FD.type = 'multi';
else
    [ix,ixc] = find(M);
    FD.fraction = [];
    FD.type = 'single';

end
FD.ix  = uint32(ix);
FD.ixc = uint32(ixc);

FD = updatetoposort(FD);
if mtf
    FD = multi_normalize(FD);
end
