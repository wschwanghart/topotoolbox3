function DEMr = resample(DEM,target,method,options)

%RESAMPLE Change spatial resolution of a GRIDobj
%
% Syntax
%
%     DEMr = resample(DEM,cellsize)
%     DEMr = resample(DEM,GRID)
%     DEMr = resample(...,method)
%
% Description
%
%     resample changes the cellsize of a grid. The function uses the
%     functions mapresize or georesize, if a target cellsize is supplied.
%     If an instance of GRIDobj is supplied as second argument, resample
%     interpolates values in DEM to match the spatial reference of GRID
%     using the function imtransform. Note that resampling might create
%     nans along the borders of the resulting raster. 
%
% Input arguments
%
%     DEM       grid object (GRIDobj)
%     cellsize  cellsize of resampled grid
%     GRID      other grid object
%     method    'cubic', 'bilinear', or 'nearest' (default is 'bilinear')
%
% Output arguments
%
%     DEMr    grid object (GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     DEMr = resample(DEM,100);
%     imagesc(DEMr)
%
% See also: griddedInterpolant, imtransform
%        
% Author:  Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 18. February, 2025 

% check also mapresize, georesize

arguments
    DEM   GRIDobj
    target  {mustBeGRIDobjOrPositiveScalar}
    method = "bilinear"
    options.crop     = true
end

method = validatestring(method,{'cubic', 'bilinear', 'nearest' },...
    'GRIDobj/resample','method',3);


% Fillvalues and methods
if isUnderlyingType(DEM,'logical') || isUnderlyingInteger(DEM)
    fillval = 0;
    method = 'nearest';
else
    fillval = nan;
end

% resample if resolution is supplied. This is easy with the 
% new functions mapresize and georesize
if ~isa(target,'GRIDobj')
    % With the new 
    if isProjected(DEM)
        scale = DEM.cellsize / target;
        [B,RB] = mapresize(DEM.Z,DEM.georef,scale,method);
    elseif isGeographic(DEM)
        scale = DEM.cellsize / target;
        [B,RB] = georesize(DEM.Z,DEM.georef,scale,method);        
    else
        if ~isempty(DEM.georef)
            scale = DEM.cellsize / target;
            [B,RB] = mapresize(DEM.Z,DEM.georef,scale,method);
        else
            % no mapping toolbox seems to be available
            error("Function requires mapping toolbox and georef information to be available.")
        end
    end
    DEMr = GRIDobj(B,RB);
else
    T = maketform('affine',[1 0 0; 0 1 0; 0 0 1]);
    switch method
        case 'cubic'
            method = 'bicubic';
    end

    if ~isEqualGeoreference(DEM,target)
        error("TopoToolbox:wronginput","Both GRIDobjs must have the same projection.")
    end

    if isProjected(DEM) || isGeographic(DEM)
        R    = DEM.georef;
        Rnew = target.georef;

        DEMr = target;
        DEMr.Z = flipud(imtransform(flipud(DEM.Z),T,method,...
        'Udata',R.XWorldLimits,'Vdata',R.YWorldLimits,...
        'Xdata',Rnew.XWorldLimits,'Ydata',Rnew.YWorldLimits,...
        'Size',Rnew.RasterSize,...
        'FillValues',fillval));
        
    else
        error("Function requires mapping toolbox")
    end

end

end
% %% ----------------------------------------------------------------------
function mustBeGRIDobjOrPositiveScalar(x)

tf = isa(x,'GRIDobj') || (isscalar(x) && x > 0);
if ~tf 
    error("TopoToolbox:validation","Input must be GRIDobj or a positive scalar.")
end
end


