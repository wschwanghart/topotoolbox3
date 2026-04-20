function DEM = pad(DEM,px,val)

%PAD Add or remove a border of pixels around a GRIDobj
%
% Syntax
%
%     DEMp = pad(DEM);
%     DEMp = pad(DEM,px,val)
%
% Description
%
%     This function adds or removes a border of constant values along all 
%     edges of a GRIDobj. By default, pad adds a one-pixel wide border with
%     zeros. To control the value and amount of padding, use the arguments
%     val and px, respectively. px can be negative. In this case, the DEM
%     is cropped by the number of pixels at each grid border.
%
%     The transformation of the grid must be 'rectilinear'.
%
% Input arguments
%
%     DEM     GRIDobj
%     px      scalar or 1x4 vector indicating the amount of padding on each 
%             of the four edges of the grid. For px=1 the resulting size of
%             DEMp is DEM.size+2. If px is a vector, then elements indicate
%             the amount of padding in [north south west east] direction.             
%     val     pad value (default=0). This value does not have any effects 
%             if px<0.
%
% Output arguments
%
%     DEMp    padded or cropped GRIDobj
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEMp20 = pad(DEM,20);
%     DEMm100 = pad(DEM,-100);
%     subplot(1,2,1)
%     imageschs(DEMp20)
%     subplot(1,2,2)
%     imageschs(DEMm100)
% 
% See also: GRIDobj/crop, padarray
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 11. June, 2024

arguments
    DEM   GRIDobj
    px    (1,4) {mustBeNumeric} = 1
    val   {mustBeScalarOrEmpty} = 0
end

% If all px are zero, we can return the original GRIDobj
if all(px == 0)
    return
end

% Make sure that we have integer values
px  = sign(px).*ceil(abs(px));

% new size of output grid
% px(1): north px(2):south px(3): west px(4): east
newsize = [DEM.size(1) + px(1) + px(2) ... 
           DEM.size(2) + px(3) + px(4)];

% check if size is negative and issue an error
if any(newsize < 0)
    error('TopoToolbox:wronginput',...
        'The number of pixels removed from the grid edges is too large')
end

% create new array
if val == 0
    Znew = zeros(newsize,class(DEM.Z));
elseif isnan(val)
    Znew = nan(newsize,class(DEM.Z));
else
    Znew = zeros(newsize,class(DEM.Z))+val;
end

% transfer values to new array
if all(px>0)
    Znew(px(1)+1:end-px(2),px(3)+1:end-px(4)) = DEM.Z;
elseif all(px<0)
    abspx = abs(px);
    Znew = DEM.Z(abspx(1)+1:end-abspx(2),abspx(3)+1:end-abspx(4));
else
    pxpos = max(px,0);
    pxneg = -min(px,0);

    Znew(pxpos(1)+1:end-pxpos(2),pxpos(3)+1:end-pxpos(4)) = ...
        DEM.Z(pxneg(1)+1:end-pxneg(2),pxneg(3)+1:end-pxneg(4));
end

% adjust referencing
DEM.Z           = Znew;
DEM.wf(:,3)     = DEM.wf(:,3) -[px(3) px(1)]'.*[DEM.wf(1,1) DEM.wf(2,2)]';
DEM.size        = size(DEM.Z);

if ~isempty(DEM.georef)
    if isProjected(DEM)
        
        R = DEM.georef;
        Rnew = maprasterref(DEM.wf,DEM.size,R.RasterInterpretation);
        Rnew.ProjectedCRS = R.ProjectedCRS;

    elseif isGeographic(DEM)
        R = DEM.georef;
        Rnew = georasterref(DEM.wf,DEM.size,R.RasterInterpretation);
        Rnew.GeographicCRS = R.GeographicCRS;

    else 
        R = DEM.georef;
        
        switch class(R) 
            case 'map.rasterref.MapCellsReference'
                Rnew = maprasterref(DEM.wf,DEM.size,R.RasterInterpretation);
            case 'map.rasterref.GeographicCellsReference'
                Rnew = georasterref(DEM.wf,DEM.size,R.RasterInterpretation);
        end
    end

    DEM.georef = Rnew;
    
    % Now check whether DEM.size == DEM.georef.RasterSize
    if ~isequal(DEM.size,DEM.georef.RasterSize)
        warning('We have a problem!')
    end
    
end

