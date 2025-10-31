function rgb = imageschs(DEM,A,options)

%IMAGESCHS Plot hillshade image with overlay
%
% Syntax
%
%     imageschs(DEM)
%     imageschs(DEM,A)
%     imageschs(_____,pn,pv,...)
%     RGB = imageschs(...)
%
% Description
%
%     Hillshading is a very powerful tool for relief depiction. imageschs
%     calculates a hillshade of a digital elevation model and colors it
%     according to a second grid or matrix. DEM must be an instance of  
%     GRIDobj and A must be a GRIDobj or matrix. 
%
%     imageschs allows to set a number of parameter name/value pairs that
%     control lighting direction and representation of missing values
%     (nan).
%
%     If called with an output variable, imageschs does not plot the
%     hillshade but returns an RGB image. The hillshading algorithm follows
%     the logarithmic approach to shaded relief representation of Katzil
%     and Doytsher (2003).
%
%     If the DEM has a geographic coordinate system (isGeographic(DEM)),
%     the function will calculate the average spatial resolution in
%     longitudinal direction and adjust DEM exaggeration accordingly. 
%
% Input
%
%     DEM         digital elevation model (GRIDobj)
%     A           coloring matrix or GRIDobj
%
%     Parameter name/value pairs
%
%     caxis or clim    two element vector defining the value range. Default 
%                      is [min(A) max(A)]. clim has precedence over caxis.  
%     colorbar         false or true (default)
%     colorbarlabel    string. Title for the colorbar
%     colorbarylabel   string. Label for the y-axis of the colorbar
%     colormap         string for colormap name or [ncol x 3] matrix. Note 
%                      that if NaNs or Infs are found in A, the colormap  
%                      must not have more than 255 colors. Default: 'jet'
%     percentclip      scalar prc (%) that truncates the displayed range to
%                      the prc's and 100%-prc's percentile of the data in
%                      A. This parameter is ignored if 'caxis' is defined.
%                      Default is prc=0. (see GRIDobj/prcclip)
%     truecolor        three element vector (rgb) with values between 0 and
%                      1 that indicates how true values are plotted if A is
%                      logical. This option also takes color names that can
%                      be interpreted by letter2rgb or ttclrr. Default is
%                      [0 1 0].
%     falsecolor       three element vector (rgb) with values between 0 and   
%                      1 that indicates how false values are plotted if A  
%                      is logical. This option also takes color names that 
%                      can be interpreted by letter2rgb or ttclrr.
%                      Default is [1 1 1].
%     nancolor         three element vector (rgb) with values between 0 and   
%                      1 that indicates how NaNs and Infs are plotted 
%                      Default is [1 1 1]. This option also takes color 
%                      names that can be interpreted by letter2rgb or 
%                      ttclrr.
%     brighten         Scalar between -1 and 1. Shift intensities of all 
%                      colors in the colormap. If <0, then colormap will be
%                      darkened, and if >0, it will be brightened.
%                      Default is 0
%     usepermanent     controls whether the hillshade is retained in 
%                      memory as persistent variable. Default is false. If 
%                      set to true, and if the DEM has the same size as the
%                      persistently stored hillshade, the function will
%                      reuse this variable, thus avoiding to recalculate
%                      hillshading.
%     medfilt          use median filter to smooth hillshading 
%                      (default=false)
%     method           'surfnorm' or 'mdow'. mdow is the multidirectional
%                      oblique hillshade algorithm.
%     azimuth          azimuth angle of illumination, (default=315)
%     altitude         altitude angle of illumination, (default=60)
%     exaggerate       elevation exaggeration (default=2). Increase to
%                      pronounce elevation differences in flat terrain. 
%     gcsadjust        {true} or false. If true, and if the DEM is in a
%                      geographic coordinate system, then exaggeration will
%                      be adjusted to the cellsize in the center of the
%                      grid.
%     ticklabels       'default', 'nice' or 'none'
%     tickstokm        true or {false}. If set to true, coordinates will be
%                      divided by 1000.
%     gridmarkers      two element vector with [dx dy] spacing of + markers
%     gridmarkercolor  three element vector (rgb) or color abbreviations
%                      as given in LineSpec (default = 'k')
%                 
% Output
%
%     RGB         [DEM.size 3] image (UINT8) with values between 0 and 255
%
%
% Example 1: Hillshade of DEM with colors derived from elevations. 
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     imageschs(DEM);
%
% Example 2: Hillshade of DEM with non-default colormap 
%     
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     imageschs(DEM,[],'colormap',landcolor);
%
% Example 3: Gray-colored hillshading only
% 
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     imageschs(DEM,[],'colormap',[.9 .9 .9],'colorbar',false);
%
% Example 4: Slope map with underlying hillshade
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     G   = gradient8(DEM);
%     imageschs(DEM,G);
%
% Example 5: Nice axis tick labels
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     imageschs(DEM,DEM,'ticklabels','nice','colorbar',true);
% 
% Example 6: Map logical GRIDobj with custom colors
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     L   = DEM>1500;
%     imageschs(DEM,L,'ticklabels','nice','colorbar',true,...
%                       'truecolor',[1 0 0],'falsecolor',[0 0 1]);
%
% Example 7: Calculate hillshade RGB and export as geotiff
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     RGB = imageschs(DEM);
%     geotiffwrite('file.tif',RGB,DEM.georef.SpatialRef,...
%            'GeoKeyDirectoryTag',DEM.georef.GeoKeyDirectoryTag);
%
% Example 8: Calculate hillshade RGB and display using imshow
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     RGB = imageschs(DEM);
%     RGB = flipud(RGB);
%     [~,R] = GRIDobj2im(DEM);
%     imshow(RGB,R)
%     axis xy
%
% References
%
%     Katzil, Y., Doytsher, Y. (2003): A logarithmic and sub-pixel approach
%     to shaded relief representation. Computers & Geosciences, 29,
%     1137-1142.
%
% See also: GRIDobj/hillshade, imagesc, ttclrr, letter2rgb,
%     GRIDobj/prcclip, adjustgeoaspectratio 
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 31. July, 2025


% Change log
%
% 23.1.2014: line 198: A = gray2ind(mat2gray(A,ncolors); 
%            replaced with
%            A = gray2ind(mat2gray(A,alims),ncolors);
% 02.6.2015: changed help and updated to R2014b
% 04.3.2016: added option medfilt and added example
% 18.3.2016: added option percentclip
% 13.6.2016: added option makepermanent
% 14.6.2016: added options 
% 06.4.2018: new example
% 04.6.2024: changed to be fit for TopoToolbox 3
% 30.8.2024: arguments block
% 10.12.2024: better handling of colors
% 31.7.2024: added bindings to libtopotoolbox


%% Input argument parsing
arguments
    DEM     GRIDobj
    A       = DEM
    options.colormap = parula(255)
    options.caxis    = []
    options.clim     = []
    options.percentclip = []
    options.truecolor   = [0 1 0]
    options.falsecolor  = [1 1 1]
    options.nancolor    = [1 1 1]
    options.exaggerate  (1,1) {mustBePositive} = 1
    options.azimuth     (1,1) {mustBeNumeric} = 315
    options.altitude    (1,1) {mustBeNumeric,mustBePositive} = 60
    options.colorbar    (1,1) = true
    options.medfilt     (1,1) = false
    options.ticklabels  = 'default'
    options.gridmarkers = []
    options.gridmarkercolor = 'k'
    options.useparallel (1,1) = true
    options.usepermanent (1,1) = false
    options.colorbarlabel = []
    options.colorbarylabel = []
    options.tickstokm (1,1) = false
    options.method = 'surfnorm'
    options.brighten {mustBeInRange(options.brighten,-1,1)} = 0
    options.gcsadjust = true
    options.uselibtt = true

end

% If the option usepermanent = true, then the function will use the
% precomputed hillshade
persistent H

% If second argument is supplied as empty array, then the DEM will be used
% to colorize the hillshade.
if isempty(A)
    A = DEM;
end

colmapfun  = options.colormap;
cbar       = options.colorbar;
truecol    = options.truecolor;
falsecol   = options.falsecolor;
exag       = options.exaggerate;
azi        = options.azimuth;
alti       = options.altitude;
nancolor   = options.nancolor;
usepermanent  = options.usepermanent;
tokm           = options.tickstokm ~= 0;
colorBarLabel  =options.colorbarlabel;
colorBarYLabel =options.colorbarylabel;
meth       = validatestring(options.method,{'default','surfnorm','mdow'});
ticklabels = validatestring(options.ticklabels,{'default','none','nice'});
gridmarkers= options.gridmarkers;
gridmarkercolor = options.gridmarkercolor;
clims      = options.caxis;
if ~isempty(options.clim)
    clims      = options.clim;
end

%% Adjust pixelsize if DEM has a geographic coordinate system
% If the DEM is in a geographic coordinate system, adjust exaggeration
% according to the pixel width in the center of the DEM.
if isGeographic(DEM) && options.gcsadjust
    [lon,lat] = getcoordinates(DEM);
    mlat = mean(lat);
    dlon = distance(mlat,lon(1),mlat,lon(2),DEM.georef.GeographicCRS.Spheroid);
    exag = exag*(DEM.cellsize)/dlon;
end
    
%% 
% check if input matrices align
validatealignment(DEM,A)

% any percentile clipping?
if ~isempty(options.percentclip)
    if ~isa(A,'GRIDobj')
        A = GRIDobj(DEM,A);
    end
    [clims,A] = prcclip(A,options.percentclip);
end

if isa(A,'GRIDobj')
    A = A.Z;
end

% constrain color range to values given in caxis
if ~isempty(clims)
    A(A<clims(1)) = clims(1);
    A(A>clims(2)) = clims(2);
end

% coordinate vector
[x,y] = getcoordinates(DEM);

% convert coordinates to km if wanted
if tokm    
    x = x*1e-3;
    y = y*1e-3;
end

% nr of colors
nhs = 256;

% calculate hillshading
if usepermanent && isequal(size(H),DEM.size)

else
    H = hillshade(DEM,'exaggerate',exag,'azimuth',azi,'altitude',alti,...
        'useparallel',options.useparallel,'method',meth,'usefused',true,...
        'uselibtt',options.uselibtt);

    H = H.Z;
    Inan = isnan(H);
    if any(Inan(:))
        H(Inan) = 1;
        clear Inan
    else
        clear Inan
    end
    
    % median filtering, if required
    if options.medfilt
        H = medfilt2(H,[3 3],'symmetric');
    end
    
    H = gray2ind(H,nhs);
end

% derive coloring
if ~isa(A,'logical')
    Inan = isnan(A(:)) | isinf(A(:));
    if any(Inan(:))
        nans = true;
        A(Inan) = nan;
    else
        nans = false;
        clear Inan
    end
    
    if ischar(colmapfun) || isstring(colmapfun)        
        ncolors = 256-nans;
        try 
            cmap = ttscm(colmapfun,ncolors);
        catch

            colmapfun = str2func(lower(colmapfun));
            cmap = colmapfun(ncolors);
        end
        
    else
        ncolors = size(colmapfun,1);
        if nans && (ncolors >= 256)
            error('TopoToolbox:GRIDobj',['NaNs found in the second input grid. \n'...
                  'Please provide colormap with less than 256 colors']);
        else
            cmap = colmapfun;
        end        
    end
    
    if cbar && isempty(clims)
        alims = [min(A(:)) max(A(:))];
    elseif cbar && ~isempty(clims)
        alims = sort(clims,'ascend');
    else
        alims = [min(A(:)) max(A(:))];
    end
    
    A = gray2ind(mat2gray(A,double(alims)),ncolors);
    
else
    ncolors = 2;
    % Convert letters to rgb
    if ischar(falsecol) || isstring(falsecol)
        falsecol = getclr(falsecol);
    end
    if ischar(truecol) || isstring(truecol)
        truecol = getclr(truecol);
    end
    cmap = [falsecol; truecol];
	nans = false;
    alims = [0 1];
end
    
% create colormap for indexing
cmap = cmap(:);
cmap = bsxfun(@times,cmap,linspace(0,1,nhs));
cmap = reshape(cmap,[ncolors 3 nhs]);
cmap = permute(cmap,[3 1 2]);
cmap = reshape(cmap,[ncolors*nhs 3]);

if options.brighten
    cmap = brighten(cmap,options.brighten);
end

% create image that indexes into the new colormap
IND  = uint16(H+1) + nhs*uint16(A) + 1;

% handle NaNs
if nans
    if ischar(nancolor) || isstring(nancolor)
        nancolor = getclr(nancolor);
    end
    cmapnan   = bsxfun(@times,nancolor,linspace(0,1,nhs)');
    IND(Inan) = uint16(H(Inan)) + nhs*(ncolors) +1;% unclear if this is ok...
    cmap      = [cmap;cmapnan];
end

% same as ind2rgb but returns a mxnx3 matrix with uint8 data
cmapUINT8 = uint8(round(cmap*256));
% see Rob Campbell's mat2im
% http://www.mathworks.de/matlabcentral/fileexchange/26322-mat2im
RGB=reshape(cmapUINT8(IND(:),:),[size(IND),3]);

% use permanent
if ~usepermanent
    H = [];
end

% plot
if nargout == 0
    imagesc(x,y,RGB);
    axis xy
    axis image
    
    % add colorbar if needed
    if cbar
        if alims(1) ~= alims(2) 
            clim(alims);
        end
        colormap(gca,cmap(nhs:nhs:nhs*ncolors,:));
        cc = colorbar;%('location','south');
        if ~isempty(colorBarLabel)
            title(cc,colorBarLabel);
        end
        if ~isempty(colorBarYLabel)
            ylabel(cc,colorBarYLabel);
        end
    end
    
    % plot nice ticklabels if wanted
    switch ticklabels
        case 'none'
            set(gca,'XTickLabel',{},'YTickLabel',{});
        case 'nice'
            
            niceticks(gca);
    end
    
    % plot grid
    if ~isempty(gridmarkers)
        if isscalar(gridmarkers)
            gridmarkers = [gridmarkers gridmarkers];
        end
          
        xgridmarkers = unique(x-rem(x,gridmarkers(1)));
        ygridmarkers = unique(y-rem(y,gridmarkers(2)));
        hold on
        [xx,yy] = meshgrid(xgridmarkers,ygridmarkers);
        plot(xx(:),yy(:),'+','Color',gridmarkercolor);
        hold off
    end
        
    
elseif nargout == 1
    rgb = RGB;
end
end


function clr = getclr(str)
% This function retrieves a color triplet based on a character or string
% input
try 
    clr = letter2rgb(str);
catch
    clr = ttclrr(str);
end
end

