function N = flowconvergence(FD)

%FLOWCONVERGENCE Compute flow convergence of a digital elevation model
% 
% Syntax
%
%     N = flowconvergence(FD)
%
% Description
%
%     flowconvergence computes the number of neighboring cells that 
%     drain into each cell. Values in N range between 0 on ridges to 8 in
%     pits. N can thus be used as a measure of local flow convergence in a
%     digital elevation model. 
%
% Input
%
%     FD      flow direction (class: FLOWobj)
%
% Output
%
%     N       grid containing the sum of contributing neighbor cells. 
%             (class: GRIDobj)
%
% Example
% 
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','c');
%     C  = flowconvergence(FD);
%     imageschs(DEM,C)
%
% See also: FLOWOBJ
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 31. August, 2024


% 4/3/2016: the function now makes copies of FD.ix and FD.ixc (see 
% FLOWobj/flowacc

ix = FD.ix;
ixc = FD.ixc;

switch FD.type
    case 'single'
        nr = zeros(FD.size,'uint8');
        for r = 1:numel(ix)
            nr(ixc(r)) = nr(ixc(r))+1;
        end
    case {'multi' 'Dinf'}
        fr = FD.fraction;
        nr = zeros(FD.size);
        for r = 1:numel(ix)
            nr(ixc(r)) = nr(ixc(r)) + fr(r);
        end
end

% Prepare output
N = GRIDobj(FD,nr);
N.name = 'flow convergence';