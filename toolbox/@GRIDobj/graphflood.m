function HWout = graphflood(DEM,P,HW,options)

%GRAPHFLOOD Graphflood algorithm as described in Gailleton et al. (2024)
%
% Syntax
%
%     HW = graphflood(DEM)
%     HW = graphflood(DEM,P,HWin,options)
%
% Description
%
%     GraphFlood 1.0 is a physics-based, efficient algorithm for
%     approximating 2D shallow‐water dynamics across digital elevation
%     models (DEMs), designed for landscape evolution and hydrological
%     modelling. It bridges the gap between simplistic drainage-area
%     methods (e.g., D8/D∞ routing) and full hydrodynamic simulations by
%     iteratively solving for flow depth and discharge on a directed
%     acyclic graph (DAG) of surface flow (Gailleton et al., 2024).
%
% Input arguments
%
%     DEM     Digital elevation model (GRIDobj)
%     P       GRIDobj or numeric scalar of precipitation rates [m/s]
%             (default is (10e-3)/3600 which is 10 mm/hour)
%     HWin    Initial water depth [m] (default is GRIDobj(DEM)) 
% 
%     Parameter name/value pairs
%     
%     'BCs'   GRIDobj with boundary conditions. Valid pixels have values of
%             1, outflow pixels have 3, and nan-pixels have 0. By default,
%             boundary pixels (incl. pixels adjacent to nan pixels) are
%             outflow pixels.
%     'dt'    time step = 1e-3 [s]. This is not simulated time as we make 
%             the steady low assumption.
%     'manning'  GRIDobj or numeric scalar of friction coefficient. 
%             Default is 0.033. 
%     'SFD'   True to compute single flow directions, False to compute 
%             multiple flow directions. Default is False.
%     'D8'    True to include diagonal paths. Default is True.
%     'N_iterations'   Number of iterations for the simulation. 
%             Default is 100.
%
% Output arguments
%
%     HW      GRIDobj with computed water depth [m]
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEM = fillsinks(DEM);
%     HW  = graphflood(DEM,(10e-3)/3600*100)
%     imageschs(DEM,HW,'colormap',flowcolor)
%
% Reference
%
%     Gailleton, B., Steer, P., Davy, P., Schwanghart, W., and Bernard, T.:
%     GraphFlood 1.0: an efficient algorithm to approximate 2D
%     hydrodynamics for landscape evolution models, Earth Surf. Dynam., 12,
%     1295–1313, https://doi.org/10.5194/esurf-12-1295-2024, 2024.
%
% See also: GRIDobj/fillsinks
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 30. July, 2025

arguments
    DEM  GRIDobj
    P    = (10e-3)/3600
    HW   = GRIDobj(DEM,'single')
    options.BCs  {validatealignment(options.BCs,DEM)} = makeBCs(DEM)
    options.dt   {mustBePositive} = 1e-3
    options.manning = 0.033
    options.SFD (1,1) = false
    options.D8  (1,1) = true
    options.N_iterations (1,1) {mustBePositive,mustBeInteger} = 100
    options.step (1,1) {mustBePositive} = single(1e-3)
end

% Elevation
DEM.Z = single(DEM.Z);

% Precipitation
if isnumeric(P) && isscalar(P)
    P = GRIDobj(DEM,'single') + single(P);
else
    validatealignment(DEM,P)
    P.Z = single(P.Z);
end 

% Water depth
if isnumeric(HW) && isscalar(HW)
    HW = GRIDobj(DEM,'single') + single(HW);
else
    validatealignment(DEM,HW)
    HW.Z = single(HW.Z);
end 

% Manning
if isnumeric(options.manning) && isscalar(options.manning)
    manning = GRIDobj(DEM,'single') + options.manning;
else
    validatealignment(DEM,options.manning)
    manning = options.manning;
end

% BCs
validatealignment(DEM,options.BCs)
BCs = options.BCs;


HWout = GRIDobj(DEM,'single');
HWout.Z = tt_graphflood(DEM.Z,HW.Z,BCs.Z,P.Z,...
    manning.Z,single(options.dt),single(DEM.cellsize),options.SFD,options.D8,...
    options.N_iterations, single(options.step));

end

function BCs = makeBCs(DEM)
% Returns a GRIDobj of boundary condition identifiers for a DEM

% Does DEM have nans
if any(isnan(DEM))
    % NaN pixels
    I = isnan(DEM);
    BCs = GRIDobj(DEM,'uint8');
    BCs.Z(I.Z) = uint8(0);
    % Boundary pixels
    I.Z = bwperim(~I.Z);
    BCs.Z(I.Z) = uint8(3);
    % Valid pixels
    I = ~isnan(DEM) & ~I;
    BCs.Z(I.Z) = uint8(1);

else
    BCs = GRIDobj(DEM,'uint8');
    BCs.Z(:,:) = uint8(3);
    BCs.Z(2:end-1,2:end-1) = uint8(1);
end
end