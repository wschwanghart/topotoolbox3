function DEM = excesstopography(DEM,maxgradient,method)

%EXCESSTOPOGRAPHY Reconstruct surface with threshold-slope surface
% 
% Syntax
%
%    DEMext = excesstopography(DEM,maxgradient)
%    DEMext = excesstopography(DEM,maxgradient,method)
%
% Description
%
%    The excess topography (Blöthe et al. 2015) is computed by solving an
%    eikonal equation (Anand et al. 2023) constrained to lie below the
%    original DEM. Where the slope of the DEM is greater than the threshold
%    slope, the eikonal solver limits the output topography to that slope,
%    but where the slope of the DEM is lower that the threshold slope, the
%    output follows the DEM.
%
%    The eikonal equation is solved using the fast sweeping method (Zhao
%    2004), which iterates over the DEM in alternating directions and
%    updates the topography according to an upwind discretization of the
%    gradient. To constrain the solution by the original DEM, the output
%    topography is initiated with the DEM and only updates lower than the
%    DEM are accepted.
%
%    The fast sweeping method ('fsm') is simpler than the fast marching
%    method ('fmm'), requires less memory, and can be faster, particularly
%    when the threshold slopes are constant or change infrequently across
%    the domain.
%
%    Excesstopography uses function in libtopotoolbox (see 
%    https://topotoolbox.github.io/libtopotoolbox/api.html). You can check
%    whether your TopoToolbox installation is linked to libtopotoolbox
%    with the function haslibtopotoolbox. If you don't have libtopotoolbox,
%    replace this function with the excesstopography function available
%    with TopoToolbox 2. 
%
% Input arguments
%
%    DEM           Digital elevation model (class: GRIDobj)
%    maxgradient   maximum gradient (scalar, GRIDobj or lithstack) provided 
%                  as radians. A lithstack object
%    method        {'fmm'} or 'fsm
%                    
% Output arguments
%
%    DEMext        reconstructed DEM. To calculate excess topography,
%                  substract the DEMext from DEM (DEM-DEMext) (class:
%                  GRIDobj)
%
% Example 1: Excesstopography with uniform gradient (20°)
%
%    DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%    DEMEXT = excesstopography(DEM,tand(20));
%    imageschs(DEM,DEM-DEMEXT)
%
% Example 2: A big pile of sand
%
%    DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%    DEM = DEM*0;
%    DEM.Z(2:end-1,2:end-1) = 3000;
%    G   = rand(DEM);
%    DEM2 = excesstopography(DEM,G);
%    surf(DEM2);
%    camlight; material dull
%
% Reference
%
%    Anand, Shashank Kumar, Matteo B. Bertagni, Arvind Singh and Amilcare
%    Porporato (2023). Eikonal equation reproduces natural landscapes with
%    threshold hillslopes. Geophysical Research Letters, 50, 21. 
%    
%    Blöthe, J.H., Korup, O., Schwanghart, W. (2015): Large landslides lie
%    low: Excess topography in the Himalaya-Karakoram ranges. Geology 43,
%    6, 523-526.
%
%    Zhao, Hongkai (2004). A fast sweeping method for eikonal equations.
%    Mathematics of Computation, 74, 250, 603-627.
%
% See also: GRIDobj/localtopography, haslibtopotoolbox
%     
% Author: Will Kearney, Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 27. January, 2025

arguments
    DEM
    maxgradient {mustBeValidGradient(maxgradient,DEM)}
    method {mustBeMember(method,{'fsm','fmm'})} = 'fmm'
end

if ~haslibtopotoolbox
    error('libtopotoolbox required to run excesstopography.')
end

% Excesstopography can be computed with lithstack to compute threshold
% angles that vary as a function of depth.
if isa(maxgradient,'lithstack')
    L = maxgradient;
    DEM.Z = tt_excesstopography_fmm3d(single(DEM.Z),...
                                      single(L.ElevStack),...
                                      single(L.Values.Sc),...
                                      single(DEM.cellsize));
    return

elseif isnumeric(maxgradient)
    if isscalar(maxgradient)
        G = repmat(single(maxgradient),size(DEM.Z));
    else
        if ~isequal(size(DEM.Z),size(maxgradient))
            error('TopoToolbox:input','Size of GRIDobj and gradient matrix are incompatible.')
        end
        G = single(maxgradient);
    end
elseif isa(maxgradient,'GRIDobj')
    G = single(maxgradient.Z);
end

switch lower(method)
    case 'fmm'
        DEM.Z = tt_excesstopography_fmm2d(single(DEM.Z),G,single(DEM.cellsize));
    otherwise
        DEM.Z = tt_excesstopography_fsm2d(single(DEM.Z),G,single(DEM.cellsize));
end
        
end


%% Input validation
function mustBeValidGradient(maxgradient,DEM)

if isa(maxgradient,'lithstack')
    L = maxgradient;
    % check whether DEM and L are compatible
    if ~isequal(size(DEM.Z),size(L.ElevStack,[1 2]))
        error('TopoToolbox:input','Size of GRIDobj and lithstack are incompatible.')
    end
    if ~ismember("Sc", L.Values.Properties.VariableNames)
        error('L.Values must contain the variable Sc.')
    end
    if any(isnan(L.Values.Sc))
        error('L.Values must not contain NaNs.')
    end

elseif isa(maxgradient,"GRIDobj")
    validatealignment(maxgradient,DEM)
    assert(~any(isnan(maxgradient)),'GRIDobj with max. gradients must not contain nans.')

elseif isnumeric(maxgradient)
    validateattributes(maxgradient,{'double','single'},...
        {"scalar","nonnan"},'excesstopography','maxgradient',2)
else
    error('Unknown input for maxgradient.')
end
end