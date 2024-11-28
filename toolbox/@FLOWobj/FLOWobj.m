classdef FLOWobj
    
%FLOWOBJ Create flow direction object
%
% Syntax
%
%     FD = FLOWobj(DEM)
%     FD = FLOWobj(DEM,'single') 
%     FD = FLOWobj(DEM,'multi' or 'dinf')
%     FD = FLOWobj(...,'pn','pv',...)
%
% Description
%
%     FLOWobj creates a flow direction object which can be used to
%     calculate terrain attributes such as flow accumulation, drainage
%     basin delineation, flow path extraction, etc. A FLOWobj is derived
%     from a digital elevation model (DEM) which is stored as a GRIDobj.
%
%     FLOWobj(DEM,'single') calculates single flow directions. This means
%     that each pixel has only one downstream neighbor. Despite these
%     limitations, 'single' will be mostly the algorithm of choice because
%     single flow directions are the basis for stream network and their
%     analysis.
%
%     FLOWobj(DEM,'multi') derives multiple flow direction (MFD) and
%     FLOWobj(DEM,'dinf') derives D infinity according to Tarboton's (1997)
%     method (Eddins 2016). Both algorithms cannot be used to calculate
%     stream networks or drainage basins, for example.
%
%     If you need to convert a transfer matrix M to a FLOWobj, then use the
%     function M2FLOWobj.
%
% Input arguments
%
%     DEM    digital elevation model (class: GRIDobj)
%     type   {'single'}, 'multi' or 'dinf'
%     
% Parameter name/value pairs   {default}
%
%     'preprocess' --  {'carve'}, 'fill', 'none'
%            set DEM preprocessing that determines flow behavior in
%            topographic depressions and flat areas 
%     'sinks' -- logical matrix same size as dem
%            true values in sinks are treated as sinks in the digital
%            elevation model and are not filled or carved, if the
%            preprocessing option fill or carve are chosen.
%     'internaldrainage' -- {false} or true
%            set this parameter value to true if flow directions should be
%            derived in the lowest, flat regions of internal drainage 
%            basins. By default, this parameter is set to false since this
%            information is usually not required and flow paths will stop
%            when entering flat, internally drained sections.
%     'cweight' --  scalar {1}
%            adjust cost matrix if preprocessing option 'carve' has been
%            chosen. 
%     'verbose' --  {false},true
%            verbose output in the command window to track computational
%            progress. Particularly interesting when working with very
%            large matrices.
%     'mex' --  {false},true (option is deprecated and will be removed
%            shortly)
%            controls if the mex routines should be called. If true, the mex
%            functions (see function compilemexfiles) must be available on
%            the search path. Using mex functions can increase the speed at
%            which an instance of FLOWobj is constructed.
%     'uselibtt' - {true},false
%            if true, then FLOWobj will be calculated using functions in
%            libtopotoolbox. If false, then functions of the image
%            processing toolbox will be used. Currently, this option only
%            applies if no sinks are to be retained.
%
% Output
%
%     FD     flow direction object (FLOWobj)
%
% Example 1
%
%     % create and evaluate flow direction object for a DEM
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM,'preprocess','carve');
%     A  = flowacc(FD);
%     imageschs(DEM,log(A));
%
% Example 2
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM,'multi');
%     A   = flowacc(FD);
%     imageschs(DEM,log(A),'colormap',flowcolor);
%
% References: 
%
%     Freeman, T. G.: Calculating catchment area with divergent flow based
%     on a regular grid, Computers & Geosciences, 17, 413–422,
%     https://doi.org/10.1016/0098-3004(91)90048-I, 1991.
%
%     Tarboton, D. G. (1997). A new method for the determination of flow 
%     directions and upslope areas in grid digital elevation models. 
%     Water Resources Research, 33(2), 309-319.
%
%     Eddins, S. (2016). Upslope area function. Mathworks File Exchange, 
%     https://www.mathworks.com/matlabcentral/fileexchange/15818-upslope-area-functions
%
%     Schwanghart, W., Groom, G., Kuhn, N. J., and Heckrath, G.: Flow
%     network derivation from a high resolution DEM in a low relief,
%     agrarian landscape, Earth Surface Processes and Landforms, 38,
%     1576–1586, https://doi.org/10.1002/esp.3452, 2013.
%
%
% See also: GRIDobj, STREAMobj, M2FLOWobj
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 3. November, 2024


properties(GetAccess = 'public', SetAccess = 'public')
    size      % size of instance of GRIDobj from which STREAMobj was derived
    type      % flow direction type (single, multi)
    ix        % [edge attribute] topologically sorted nodes (givers)
    ixc       % [edge attribute] topologically sorted nodes (receivers)
    fraction  % [edge attribute] fraction transfered between nodes
    cellsize  % cellsize of the grid (scalar)
    wf        % 2-by-3 affine transformation matrix (see worldFileMatrix)
    georef    % additional information on spatial referencing
    
end

properties(GetAccess = 'public', SetAccess = 'public')
    fastindexing = false; % set to true to initiate fast indexing
    ixcix     % indexing matrix for fast indexing
end

methods
    function FD = FLOWobj(DEM,type,options)

        arguments
            DEM = []
            type {mustBeTextScalar} = 'single'
            options.preprocess {mustBeTextScalar,...
                mustBeMember(options.preprocess,{'carve','fill','none'})} ...
                                    = 'carve'
            options.verbose (1,1) = false
            options.sinks = []
            options.internaldrainage (1,1) = false
            options.uselibtt (1,1) = true && haslibtopotoolbox
            options.tweight (1,1) {mustBeNumeric,mustBePositive} = 2
            options.cweight (1,1) {mustBeNumeric,mustBePositive} = 1
            options.mex (1,1) = false
        end

        if nargin == 0
            % Create an empty FLOWobj
            return
        end
        
        %% Some additional input checking
        % If there are input arguments, the first input must be a GRIDobj
        assert(isa(DEM,'GRIDobj'),'TopoToolbox:WrongInput', ...
            'First input argument must be a GRIDobj')
        type = validatestring(type,{'single','multi','dinf'},'FLOWobj','type',2);


        % Transfer DEM (GRIDobj) properties to FLOWobj
        FD.cellsize = DEM.cellsize;
        FD.wf       = DEM.wf;
        FD.georef   = DEM.georef;
        FD.size     = DEM.size;
        FD.type     = type;
        nrc         = numel(DEM.Z);

        % Compute auxiliary topography
        [D,DEM,SILLS] = createAuxiliaryTopo(DEM,...
            'internaldrainage',options.internaldrainage,...
            'preprocess',options.preprocess,...
            'verbose',options.verbose,...
            'uselibtt',options.uselibtt);

        % D     - Auxiliary elevations starting at presillpixels. All other
        %         pixels are -inf (numeric matrix)
        % DEM   - Filled DEM (GRIDobj)
        % SILLS = logical matrix with locations of sills
        
        switch type
            case 'single'  
                  
                if options.mex
                    % mexed version: 
                    % identifies steepest downward neighbors in the DEM and 
                    % the distance grid obtained from graydist and performs
                    % a topological sort                   
                    D(isinf(D)) = inf;
                    D(SILLS)    = 0;
                    clear SILLS
                    dem = double(DEM.Z);
                    % call to steepest neighbor. 
                    SE = steepestneighbor_mex(dem,D,FD.cellsize);
                    if options.verbose
                       disp([char(datetime("now"))  ' -- Steepest neighbor identified'])
                    end
                    clear D
                    I  = isnan(dem);
                    clear dem
                    IX = find(SE==0 & ~I);
                    MV = nnz(I);
                    [FD.ix,FD.ixc] = tsort_mex(SE,uint32(IX),uint32(MV));
                    I  = FD.ix==0;
                    FD.ix(I) = [];
                    FD.ixc(I) = [];
                    return
                    
                else
                    %% New version (uncomment to run with this code) 
                    % But does not yet work properly because it creates
                    % some weird chessboard patterns in drainage basins and
                    % is generally not consistent with the legacy method.
                    % Also not tested fully yet, it also seems slower.
                    % -- uncomment from here to test.
                    % D(SILLS)    = 0;
                    % [ix,ixc]=ixdsneighbors_steepest(DEM.Z,D,DEM.cellsize);
                    % 
                    % FD.ix = uint32(ix);
                    % FD.ixc = uint32(ixc);
                    % FD.fraction = [];
                    % FD = updatetoposort(FD);
                    % return
                    % -- uncomment til here

                    %% new version until here

                    % Find steepest neighbor
                    % sort pixels in D
                    % adapted from sortrows.m
                    [~,IXSortedFlats] = sort(D(:),'descend');
                    clear D

                    if options.verbose
                        disp([char(datetime("now"))  ' -- Pixels sorted (1)'])
                    end

                    ndx = (uint32(1):uint32(nrc))';
                    ndx = ndx(IXSortedFlats);
                    clear IXSortedFlats

                    [~,FD.ix] = sort(DEM.Z(ndx),'descend');

                    if options.verbose
                        disp([char(datetime("now"))  ' -- Pixels sorted (2)'])
                    end

                    FD.ix = uint32(FD.ix);
                    FD.ix = ndx(FD.ix);
                    clear ndx;

                    % a fast solution that has quite much memory overhead...
                    pp = zeros(FD.size,'uint32');
                    IX = (uint32(1):uint32(numel(DEM.Z)))';
                    pp(FD.ix) = IX;

                    % cardinal neighbors
                    IXC1 = imdilate(pp,[0 1 0; 1 1 1; 0 1 0]>0);
                    xxx1 = IXC1;
                    IX   = IXC1(FD.ix);
                    IXC1 = FD.ix(IX);
                    G1   = (DEM.Z(FD.ix)-DEM.Z(IXC1))/(FD.cellsize);
                    G1(FD.ix == IXC1) = -inf;

                    % diagonal neighbors
                    IXC2 = imdilate(pp,[1 0 1; 0 1 0; 1 0 1]>0);
                    xxx2 = IXC2;
                    IX   = IXC2(FD.ix);
                    IXC2 = FD.ix(IX);
                    G2   = (DEM.Z(FD.ix)-DEM.Z(IXC2))/(norm([FD.cellsize,FD.cellsize]));

                    % choose the steeper one
                    I    = G1<=G2 & xxx2(FD.ix)>xxx1(FD.ix);
                    FD.ixc = IXC1;
                    FD.ixc(I) = IXC2(I);

                    I = FD.ixc == FD.ix;
                    FD.ix(I) = [];
                    FD.ixc(I) = [];

                    % remove nans
                    I = isnan(DEM.Z);
                    FD.ixc(I(FD.ix)) = [];
                    FD.ix(I(FD.ix)) = [];

                end

                if options.verbose
                    disp([char(datetime("now"))  ' -- Ordered topology established'])
                end
                
            
        case 'multi'
            %% Multiple flow direction
            D(SILLS) = 0;
            
            % This function returns a list of downstream neighbors ixc for 
            % each pixel ix. The list is not topologically sorted and the
            % fractions are still to be de
            [ix,ixc]=ixdsneighbors(DEM.Z,D,'keepequal',false);

            FD.ix = ix;
            FD.ixc = ixc;
            FD.fraction = double((DEM.Z(ix)-DEM.Z(ix))./ ...
                getdistance(ix,ixc,DEM.size,DEM.cellsize,'single'));
            FD.fraction(FD.fraction == 0) = 1;
            FD = updatetoposort(FD);
            FD = multi_normalize(FD);
            return
      
        case 'dinf'
            % This part may need some additional work
            R = dem_flow(DEM.Z);
            M = flow_matrix(DEM.Z,R,DEM.cellsize,DEM.cellsize);
            M = -M';

            if ~isempty(DEM.georef)
                FD = M2FLOWobj(M,DEM.georef);
            else
                FD = M2FLOWobj(M,DEM.wf);
            end
            FD.type = 'Dinf';
            
        end 
        
    end

    function FD = saveobj(FD)
        FD.fastindexing = false;
    end
    
    function FD = set.fastindexing(FD,val)
        % fastindexing enables quick traversal along flow directions
        % starting from a single pixel. It is used when deriving flow 
        % paths, e.g., flowpathextract.
        
        % switch lower(FD.type)
        %    case {'multi','dinf'}
        %        error('Fast indexing only possible for single flow directions.')
        %end

        validateattributes(val,{'numeric','logical'},{'scalar'})
        
        if val
            if ~strcmp(FD.type,'single')
                error('TopoToolbox:fastindexing','Fast indexing is only possible for single flow directions');
            end
            FD.fastindexing = true;
            FD.ixcix  = zeros(FD.size,'uint32');
            FD.ixcix(FD.ix) = uint32(1):uint32(numel(FD.ix));
        else
            FD.fastindexing = false;
            FD.ixcix  = [];
        end
    end
    
    function FD = multi_normalize(FD)
        
        if isempty(FD.fraction)
            return
        end
        s = accumarray(FD.ix,FD.fraction,[prod(FD.size) 1],@sum);
        FD.fraction = FD.fraction./s(FD.ix);
    end
             
end
end


