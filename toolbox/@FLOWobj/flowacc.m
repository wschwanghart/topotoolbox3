function OUT = flowacc(FD,W0,RR,options)

%FLOWACC flow accumulation (upslope area, contributing area)
%
% Syntax
%
%     A = flowacc(FD)
%     A = flowacc(FD,W0)
%     A = flowacc(FD,W0,RR)
%
% Description
%
%     flowacc calculates the number of upstream cells based on a flow
%     direction object (FLOWobj). To obtain upslope or contributing area in
%     metric units, multiply A by the square of the cellsize, e.g.
%     A = flowacc(FD).*(FD.cellsize^2).
%
%     The second input argument can be used to define spatially variable
%     weights into the flow accumulation e.g., to simulate spatially
%     variable precipitation patterns. By default, W0 is a grid of ones.
%
%     The third input argument is the runoff ratio. By default, the runoff
%     ratio equals one everywhere. To simulate infiltration or channel
%     transmission losses, values between 0 and 1 indicate the proportion
%     of flow that is transferred along a distance of 1 map unit. The 
%     amount of water (or other matter) transferred from one cell ix to its
%     downstream neighbor ixc is then calculated as 
% 
%     A(ixc)_1 = A(ixc)_0 + fraction(ix_ixc) * A(ix)*exp(-(1-RR(ix))*dx);
%
%     where fraction is one for single flow directions or between 0 and 1
%     for multiple flow directions, RR is the runoff ratio, and dx is the
%     distance between the pixel ix and ixc. 
%
% Input arguments
%
%     FD    Flow direction object (class: FLOWobj)
%     W0    weight grid (class: GRIDobj) 
%     RR    runoff ratio grid (class: GRIDobj) or scalar
%
% Output arguments
%
%     A     flow accumulation grid (class: GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     A = flowacc(FD);
%     imageschs(DEM,dilate(sqrt(A),ones(5)),'colormap','flowcolor')
%     
% 
% See also: FLOWobj, GRIDobj, FLOWobj/drainagebasins
% 
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
%         William Kearney
% Date: 25. June, 2026
  
arguments
    FD   FLOWobj
    W0   = []
    RR   = []
    options.uselibtt (1,1) = true
end

if ~isempty(W0)
    validatealignment(FD,W0);
end

if ~isempty(RR)
    if isa(RR,'GRIDobj')
        validatealignment(FD,RR);
    end
end

if options.uselibtt && haslibtopotoolbox
    if isempty(W0)
        W = ones(FD.size);
    elseif isa(W0, "GRIDobj")
        W = double(W0.Z);
    else
        W = double(W0);
    end
    switch (lower(FD.type))
        case 'single'
            R = ones(numel(FD.ix), 1);
        case {'multi', 'dinf'}
            R = double(FD.fraction);
    end

    R = calcEdgeTransfer(R,RR,FD);

    A = tt_traverse_down_f64_add_mul(W, R, int64(FD.ix - 1), int64(FD.ixc-1));
else
    if nargin == 1 || (nargin > 1 && isempty(W0))
        A = ones(FD.size);
    else
        if isa(W0,'GRIDobj')
            A = double(W0.Z);
        else
            A = double(W0);
        end
    end

    % copies of ix and ixc to increase speed with 2015b
    ix = FD.ix;
    ixc = FD.ixc;

    switch lower(FD.type)
        case 'single'
            if nargin < 3

                for r = 1:numel(ix)
                    A(ixc(r)) = A(ix(r))+A(ixc(r));
                end
            else

                R = calcEdgeTransfer(1,RR,FD);

                for r = 1:numel(ix)
                    A(ixc(r)) = A(ixc(r)) + R(r)*A(ix(r));
                end

            end

        case {'multi','dinf'}
            fraction = FD.fraction;
            if nargin < 3
                for r = 1:numel(ix)
                    A(ixc(r)) = A(ix(r))*fraction(r) + A(ixc(r));
                end
            else
                R = calcEdgeTransfer(FD.fraction,RR,FD);
                for r = 1:numel(ix)
                    A(ixc(r)) = R(r).*A(ix(r)) + A(ixc(r));
                end
            end
    end
end


%% Prepare Output
OUT = GRIDobj(FD,A);
OUT.zunit = 'nr of cells';
OUT.name  = 'flow accumulation';



end

function R = calcEdgeTransfer(R,RR,FD)

if isempty(RR)
    % Do nothing
elseif isscalar(RR) && isnumeric(RR)
    R = R .* exp(-(1-RR).*getdistance(FD.ix,FD.ixc,FD.size,FD.cellsize));
    R = double(R);
else 
    if isa(RR, "GRIDobj") 
        validatealignment(FD,RR)
        RR = RR.Z;
    else
        if ~isequal(FD.size,size(RR))
            error(['RR must be a matrix with ' FD.size(1) ' rows and ' FD.size(2) ' columns.'])
        end
    end
    R = R .* exp(-(1-RR(FD.ix)).*getdistance(FD.ix,FD.ixc,FD.size,FD.cellsize));
    R = double(R);
end

end
