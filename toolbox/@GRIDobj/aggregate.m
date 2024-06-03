function C = aggregate(DEMlowres,DEMhighres,options)

%AGGREGATE Resampling a GRIDobj using aggregation/binning
%
% Syntax
%
%     C = aggregate(B,A)
%     C = aggregate(B,xy)
%     C = aggregate(B,xyz)
%     C = aggregate(...,pn,pv,...)
%
% Description
%
%     This function resamples the grid A to C to match the extent and
%     resolution of grid B. B must spatially overlap with A and must have a
%     coarser resolution. By default, aggregate uses the mean to calculate
%     new values in grid C, but aggregate takes any other function that takes 
%     a vector and returns a scalar (e.g. median, std, ...).
%
%     Values to be aggregated can also be supplied as list of coordinates
%     (and attributes). This is particularly useful if point density for
%     each pixel is greater than one. This binning technique may not be 
%     appropriate if the point density is low, because this will produce a
%     grid which is only sparsely populated with values.  
%
% Input arguments
%
%     B         low resolution GRIDobj. A and B must have the same 
%               coordinate system
%     A         high resolution GRIDobj
%     xy        two column matrix with coordinates. In this case, 
%               aggfun = @numel and fillval = 0.
%     xyz       three column matrix with coordinates and third column being
%               some attribute (e.g. elevation, ...)
%
%     Parameter name/value pairs
% 
%     aggfun    anonymous function that aggregate values. The function must
%               take a vector and return a scalar. The default is @mean.
%               Other possible functions are @numel to obtain counts,
%               @median, @std, ...
%
%     fillval   pixel values in the output grid with no overlapping values 
%
% Output arguments
%
%     C         GRIDobj with same extent and resolution as B and aggregated
%               values derived from A
%
%
% See also: GRIDobj/resample, GRIDobj/reclabel, accumarray
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 3. June, 2024

arguments
    DEMlowres GRIDobj
    DEMhighres 
    options.aggfun = @mean
    options.fillval = nan
    options.class   = underlyingType(DEMhighres)
end

aggfun  = options.aggfun;
fillval = options.fillval; 
outclass   = options.class;

if isa(DEMhighres,'GRIDobj')
    % A is a high res GRIDobj
    [x,y] = getcoordinates(DEMhighres);
    IX    = zeros(DEMhighres.size);
    [X,Y] = getcoordinates(DEMlowres);

    sy = size(y);
    for r = 1:numel(x)
        IX(:,r) = coord2ind(X,Y,repmat(x(r),sy),y);
    end
    
    Z = DEMhighres.Z(:);
    IX = IX(:);
    inan = isnan(IX);
    IX(inan) = [];
    Z(inan) = [];
    
else
    % A is a list of coordinates and attributes
    xyz = DEMhighres;
    IX  = coord2ind(DEMlowres,xyz(:,1),xyz(:,2));
    inan = isnan(IX);
    xyz(inan,:) = [];
    IX(inan) = [];
    
    if size(xyz,2) == 3
        Z = xyz(:,end);
    else
        Z = true; %true(size(IX));
        aggfun = @numel;
        fillval = 0;
    end
end
    
if isempty(Z)
    C = GRIDobj(DEMlowres,outclass);
    C.Z(:,:) = cast(fillval,outclass);
    return
end

z = accumarray(IX,Z,[prod(DEMlowres.size) 1],...
    @(x) cast(aggfun(x),outclass),cast(fillval,outclass));
C = DEMlowres;
C.Z = reshape(z,DEMlowres.size);


