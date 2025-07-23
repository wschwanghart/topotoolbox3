function [D,L] = distance(DEM,x,y)

%DISTANCE distance transform
%
% Syntax
%
%     D = distance(DEM,ix)
%     D = distance(DEM,I)
%     D = distance(DEM,S)
%     D = distance(DEM,x,y)
%     D = distance(DEM,MS)
%     [D,L] = ...
%
% Description
%
%     This function calculates the distance of all pixels in the GRIDobj 
%     DEM to the nearest location in ix (linear index to pixels in DEM), I
%     (logical GRIDobj aligned with DEM), S (STREAMobj), x and y
%     coordinates, MS (mapping structure array).
%
%     Note that this function uses the Image Processing Toolbox function
%     bwdist. If near locations are provided as x-y coordinates or mapping
%     structure, the function snaps coordinates to the nearest pixels and
%     calculates distances from these pixels. 
%
% Input arguments
%
%     DEM      GRIDobj
%     ix       linear index into DEM
%     I        logical GRIDobj spatially aligned with DEM
%     S        STREAMobj spatially aligned with DEM
%     x,y      coordinates vectors
%     MS       point or line mapping structure array (as obtained from the
%              function shaperead)
%
% Output arguments
%   
%     D        distance in map units
%     L        linear index of nearest location
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'preprocess','carve');
%     S   = STREAMobj(FD,'minarea',1000);
%     D = distance(DEM,S);
%     imageschs(DEM,D)
%     hold on
%     plot(S,'w');
%     hold off
%
% See also: FLOWobj/flowdistance, bwdist
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 11. July, 2025

arguments
    DEM GRIDobj
    x = []
    y = []
end

if nargin == 1
    ix = find(DEM);
elseif nargin == 2
    if isnumeric(x)
        ix = x;
    elseif isa(x,'GRIDobj')
        ix = x.Z~=0;
    elseif isa(x,'STREAMobj')
        ix = x.IXgrid;
    elseif isstruct(x)
        MS = x;
        
        switch lower(MS(1).Geometry)
            case 'point'
                % possible fieldnames
                possfnx = {'x' 'X' 'lon' 'Lon'};
                possfny = {'y' 'Y' 'lat' 'Lat'};
                for counter = 1:numel(possfnx)
                    
                    try
                        x = [MS.(possfnx{counter})];
                        y = [MS.(possfny{counter})];
                        
                        break
                        
                    catch
                        continue
                    end
                
                end
                ix = coord2ind(DEM,x,y);
            case {'line' 'polyline'}
                ix  = line2GRIDobj(DEM,x);
                ix  = ix.Z;
            otherwise 
                error('unsupported format of the structure array')
        end
    else
        error('second input has unsupported format');
    end
    
                
else
    if ~isequal(size(x),size(y))
        error('The second and third input must have the same size.')
    end
    ix = coord2ind(DEM,x,y);
end

MASK = false(DEM.size);
MASK(ix) = true;

D = GRIDobj(DEM);

if nargout == 1
    D.Z = bwdist(MASK,'e');
else    
    [D.Z,L] = bwdist(MASK,'e');
end
D = D.*DEM.cellsize;
                
        
        
        
    
    
