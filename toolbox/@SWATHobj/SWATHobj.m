classdef SWATHobj
%SWATHOBJ Create swath profile object (SWATHobj)
%
% Syntax
%
%    SW = SWATHobj(DEM)
%    SW = SWATHobj(DEM,x,y)
%    SW = SWATHobj(DEM,x,y,d)
%    SW = SWATHobj(...,'pn','pv')
%
%
% Description
%
%     SWATHobj creates a swath profile object, which can be used to
%     obtain statistics on terrain attributes such as elevation or slope
%     angles along a line or other directional x,y pairs.
%
%     SWATHobj(DEM) opens a figure to interactively create a SWATHobj with
%     the user defining a line with an arbitrary amount of nodes.
%
%     SWATHobj(DEM,x,y) creates a SWATHobj from x,y vectors.
%
%     SWATHobj(DEM,x,y,d) creates a SWATHobj from x,y vectors and
%     the distance vector d, which has to be of the same size as x,y.
%
%     SWATHobj(DEM,'pn','pv'...) creates a SWATHobj and defines certain
%     parameter values that control the geometry of the SWATHobj.
%
%
% Input arguments
%
%     DEM    digital elevation model (Class: GRIDobj)
%     xy     directional x,y data pair (n x 2 array)
%
% Parameter name/value pairs   {default}
%
%     'width'    scalar {1e4}
%            width of the swath profile in meters
%
%     'gap'    scalar {0}
%            width of a gap centered on the trace of the swath profile
%            within which no data is obtained
%
%     'dx'    scalar {cellsize of DEM}
%            resampling distance in the longitudinal direction of the swath
%            profile. Provide empty matrix ([]) for no resampling at all.
%
%     'dy'    scalar {cellsize of DEM}
%            resampling distance in the transverse direction of the swath
%            profile. Provide empty matrix ([]) for no resampling at all.
%
%     'keepnodes'    {false},true
%            switch to adjust resampling of points along swath profile to
%            make sure the original nodes are included. If activated (true)
%            spacing of points along profile will most likely not be
%            unique.
%
%     'keepdist'    false,{true}
%            switch to adjust distance vector of swath profile to match
%            original distance vector of input data. If deactivated (false)
%            distances will change according to resampling and smoothing.
%
%     'keeptrace'    {false},true
%            switch to adjust trace of swath profile to match original
%            trace of input data. If activated (true) but the trace of the
%            profile has been resampled and/or smoothed, ...
%
%     'smooth'    {false}, true
%            optional smoothing of profile trace in x,y space. The function
%            uses smoothdata to smooth the profile trace
%
%     'smoothmethod' {'sgolay'}
%            method to smooth the profile trace. See smoothdata for more
%            options.
%
%     'smoothingfactor' value between 0 and 1 (default = 0.25)
%            see smoothdata for explanation of the smoothing factor
%
%     'hillshade' {false},true
%            use hillshade when interactively drawing the profile trace.
%
% Output
%
%     SW     swath profile object (SWATHobj)
%
%
% SWATHobj properties:
%
%     xy0       - raw x,y points used to create SWATHobj
%     zd0       - elevation and distance of raw points
%     dx        - longitudinal resampling interval (xy-unit of DEM)
%     dy        - transverse resampling interval (xy-unit of DEM)
%     width     - width of SWATHobj (xy-unit of DEM)
%     gap       - central gap along SWATHobj (xy-unit of DEM)
%     smooth    - cell with scalars corresponding to smoothing values
%     xy        - x,y points of central line of SWATHobj
%     distx     - distance along SWATHobj
%     disty     - distance across SWATHobj
%     X         - X coorindates of SWATHobj data points
%     Y         - Y coorindates of SWATHobj data points
%     Z         - Z values of SWATHobj at data points
%     name      - name of SWATHobj
%     xyunit    - xyunit taken from DEM (GRIDobj)
%     zunit     - zunit taken from DEM (GRIDobj)
%     georef    - georeference structure taken from DEM (GRIDobj)
%
%
% Example 1
%
%     % interactively create swath profile
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     SW = SWATHobj(DEM)
%
% Example 2
%
%     % create SWATHobj along a STREAMobj
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     A  = flowacc(FD);
%     S = STREAMobj(FD,A>100);
%     S = trunk(klargestconncomps(S,1));
%     [x,y] = STREAMobj2XY(S);
%     ix = ~isnan(x);
%     SW = SWATHobj(DEM,x(ix),y(ix),'smooth',200);
%
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de)
%         Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: May, 2024
    
    
    properties
        xy0     % add descriptions
        zd0
        dx
        dy
        width
        gap
        smooth
        xy
        distx
        disty
        X
        Y
        Z
        name
        xyunit
        zunit
        georef
    end
    
    methods
        
        function SW = SWATHobj(DEM,x,y,d,options)

            arguments
                DEM   GRIDobj
                x = []
                y = []
                d = []
                options.width (1,1) = 1e4
                options.dx (1,1) = DEM.cellsize
                options.dy (1,1) = DEM.cellsize
                options.gap (1,1) = 0
                options.keepdist (1,1) = true
                options.keeptrace (1,1) = false
                options.keepnodes (1,1) = false
                options.smooth (1,1) = false
                options.smoothmethod = 'sgolay'
                options.hillshade (1,1) = false
                options.smoothingfactor (1,1) = .25
            end

            if isa(x,'STREAMobj')
                error(['STREAMobj as input in standard call to SWATHobj ' ...
                        'not supported. Use STREAMobj2SWATHobj instead.'])
            end
            
            drawtrace = false;
            if isempty(x) % SWATHobj is created interactively
                drawtrace = true;
                hfig = figure;
                ax   = axes(hfig);
                if options.hillshade
                    imageschs(DEM,[],'colorbar',false)
                else
                    imagesc(DEM,'parent',ax)
                end
                title('Draw a line (double-click to end)')
                [XY] = getline;
                x0 = XY(:,1);
                y0 = XY(:,2);
                d0 = getdistance(x0,y0);
                hold on; 
                plot(x0,y0,'k-'); 
                plot(x0,y0,'w--');
                hold off
            elseif ~isempty(x) && isempty(y)
                validateattributes(x,{'numeric'},{'2d','ncols',2})
                x0 = x(:,1);
                y0 = x(:,2);
                d0 = getdistance(x0,y0);
            elseif ~isempty(x) && ~isempty(y)
                validateattributes(x,{'numeric'},{'vector'});
                validateattributes(y,{'numeric'},{'vector','size',size(x)});
                x0 = x(:);
                y0 = y(:);
                if isempty(d)
                    d0 = getdistance(x0,y0);
                else
                    validateattributes(d,{'numeric'},{'vector','size',size(x)});
                    if ~all(diff(d)>0) || ~all(diff(d)<0)
                        error('Distances (3rd input variable) must be increasing or decreasing')
                    end
                    d0 = d;
                end
            end

            % construct SWATHobj
            SW.width    = options.width;
            SW.gap      = options.gap;
            SW.dx       = options.dx;
            SW.dy       = options.dy;
            SW.georef   = DEM.georef;
            SW.name     = ['SWATHobj created from GRIDobj ',DEM.name];
            SW.zunit    = DEM.zunit;
            SW.xyunit   = DEM.xyunit;
            
            keepdist    = options.keepdist;
            keepnodes   = options.keepnodes;
            keeptrace   = options.keeptrace;
            smoothing   = options.smooth;
            
            
            %% Create SWATHobj
            if ~isempty(SW.dx) % resample along x
                if keepnodes % SIMPLIFY??
                    this_x = []; this_y = []; this_d = [];
                    for k = 1 : length(x0)-1
                        [xt,yt,dt] = interpline(x0(k:k+1),y0(k:k+1),d0(k:k+1),SW.dx);
                        this_x = [this_x;xt];
                        this_y = [this_y;yt];
                        this_d = [this_d;d0(k)+dt];
                    end
                else
                    [this_x,this_y,this_d] = interpline(x0,y0,d0,SW.dx);
                end
            else % no resampling along x
                this_x = x0;
                this_y = y0;
            end
            ix0 = coord2ind(DEM,x0,y0);
            z0 = DEM.Z(ix0);
            
            if length(this_x)>1
                SW.xy0 = [x0(:),y0(:)];
                SW.zd0 = [z0(:),d0(:)];

                if smoothing>0 % use moving average?

                    temp = smoothdata([this_x(:) this_y(:)],...
                        options.smoothmethod,...
                        'SmoothingFactor',options.smoothingfactor);
                    this_xf = temp(:,1);
                    this_yf = temp(:,2);
                    SW.smooth = options.smoothingfactor;

                else
                    this_xf = this_x;
                    this_yf = this_y;
                end

                dX = diff(this_xf); % dx between points along profile
                dY = diff(this_yf); % dy between points along profile
                
                if ~keepdist
                    this_d = getdistance(this_xf,this_yf);
                end
                
                if ~keeptrace
                    this_x = this_xf;
                    this_y = this_yf;
                end
                
                SW.xy = [this_x,this_y];
                SW.distx = this_d;
                
                hwidth = SW.width/2;
                hgap = SW.gap/2;
                if isempty(SW.dy)
                    SW.disty = -hwidth;
                else
                    SW.disty = (-hwidth : SW.dy : -hgap)';
                end
                SW.disty = [SW.disty; flipud(abs(SW.disty))];
                SW.disty = unique(SW.disty,'stable');
                
                % Create transverse profiles that are orthogonal to the center line
                DX = ([dX(1); dX]+[dX; dX(end)])./2;
                DY = ([dY(1); dY]+[dY; dY(end)])./2;
                [theta,~] = cart2pol(DX,DY);  % direction between points along profile
                theta_orthogonal = theta+pi/2; % orthogonals to center line of swath profile
                [x_orthogonal,y_orthogonal] = pol2cart(theta_orthogonal,1); % dx, dy of orthogonals
                % create new points
                ny = length(SW.disty);
                nx = length(SW.distx);
                SW.X = repmat(SW.disty,1,nx) .* repmat(x_orthogonal',ny,1) + repmat(this_x',ny,1);
                SW.Y = repmat(SW.disty,1,nx) .* repmat(y_orthogonal',ny,1) + repmat(this_y',ny,1);
                
                % Interpolate GRIDobj values
                SW.Z = nan(size(SW.X));
                ix = 1:numel(SW.X);
                SW.Z(ix) = interp(DEM,SW.X(ix),SW.Y(ix));
            else
                error('Too few points to create SWATHobj.')
            end
            
            if drawtrace
                hold on
                plot(this_xf,this_yf,'b-')
                plot(this_xf,this_yf,'w--')
                hold off
            end
        end
        
    end % methods
    
end % classdef
