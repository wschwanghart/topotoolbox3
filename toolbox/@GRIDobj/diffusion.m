function DEM = diffusion(DEM,options)

%DIFFUSION Solve the diffusion equation
%
% Syntax
%
%     DEMd = diffusion(DEM)
%     DEMd = diffusion(DEM,pn,pv,...)
%
% Description
%
%     This function uses an implicit scheme to solve the linear diffusion
%     equation for a DEM.
%
% Input arguments
%
%     DEM     digital elevation model
%
%     Parameter name/value pairs
%
%     D           diffusivity (m^2 /y) (default = 1)
%     timespan    duration (y) (default = 1000)
%     numsteps    number of iterations (default = 5)
%     streamnet   STREAMobj
%     solver      'pcg' (default) or '\'
%     pcgtol      1e-6 (default)
%     uplift      GRIDobj or scalar (default = 0) in mm/y
%
% Output arguments
%
%     DEMd    diffused DEM
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEMd = DEM;
%     for r = 1:10;
%     DEMd = diffusion(DEMd); 
%     imageschs(DEMd); 
%     drawnow; 
%     end
%     figure
%     imageschs(DEM,DEMd-DEM)
%
% See also: GRIDobj/filter, ttlem
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 12. June, 2026

arguments
    DEM
    options.timespan (1,1) {mustBeNumeric,mustBePositive} = 1000
    options.numsteps (1,1) {mustBeInRange,mustBePositive} = 5
    options.D = 1
    options.streamnet = []
    options.solve {mustBeMember(options.solve,{'\','pcg'})} = 'pcg'
    options.pcgtol (1,1) {mustBeNumeric,mustBePositive} = 1e-6
    options.uplift = 0
end

D = options.D;
numsteps = options.numsteps;
dt = options.timespan / numsteps;

nrc = prod(DEM.size);

% get neighbor indices
[ic,icd] = ixneighbors(DEM.Z,[],4);

% remove indices to nan-cells
I = isnan(DEM.Z(ic)) | isnan(DEM.Z(icd));
ic(I) = [];
icd(I) = [];

% calculate laplacian
L        = sparse(ic,icd,1,nrc,nrc);
L        = spdiags(sum(L,2),0,nrc,nrc) - L;

% and the diffusion matrix
if isa(D,'GRIDobj')
	D = D.Z(:);
else
	D = repmat(D,numel(DEM.Z),1);
end

if isempty(options.streamnet)
	D  = speye(nrc) + spdiags(D,0,nrc,nrc)*dt/(2*DEM.cellsize^2)*L;
else
	S  = +STREAMobj2GRIDobj(options.streamnet);
	D  = speye(nrc) + spdiags(D,0,nrc,nrc)*dt/(2*DEM.cellsize^2)*spdiags(1-S.Z(:),0,nrc,nrc)*L;
end

% must work with doubles
% remember class
c        = class(DEM.Z);
% ensure double
DEM.Z    = double(DEM.Z);

% set nan values to zero
I        = isnan(DEM.Z);
DEM.Z(I) = 0;

% solve
Z1 = DEM.Z(:);

if isa(options.uplift,'GRIDobj')
    u = options.uplift.Z;
else 
    u = double(options.uplift);
    u = repmat(u,DEM.size);
end

if ~isempty(options.streamnet)
    u(S.Z>0) = 0;
end
u = u(:);
u = u/1000 * dt;

for r = 1:numsteps
    switch options.solver
        case '\'
            Z1 = D\(Z1+u);
        case 'pcg'
            [Z1,~] = pcg(D,double(Z1+u),options.pcgtol,[],[],[],Z1);
        otherwise
            error('unknown solver')
    end
end

% reshape
DEM.Z = reshape(Z1,DEM.size);
% reset values to nan
DEM.Z(I) = nan;
% reset input class
DEM.Z = cast(DEM.Z,c);