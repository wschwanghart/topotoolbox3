function OUT = mapswath(SW,M,method)
%MAPSWATH Obtains Z values along a SWATHobj from an arbitrary GRIDobj
%
% Syntax
%
%     OUT = mapswath(SW,G)
%     OUT = mapswath(SW,G,method)
%     OUT = mapswath(SW,S,nal)
%
% Description
%
%     MAPSWATH interpolates Z values from the GRIDobj M, at the data point
%     locations of the SWATHobj SW. If SW and M are of different spatial
%     reference (projection), MAPSWATH tries to use the mapping toolbox to 
%     convert between units.
%
% Input arguments
%
%     SW      instance of SWATHobj
%     G       instance of GRIDobj
%     method  interpolation method, can be any of 'linear','nearest',
%             'spline','pchip','cubic'; by default 'linear'
%     S       instance of STREAMobj
%     nal     node-attribute list
%
% Output arguments
%
%     OUT    instance of SWATHobj
%
% Examples
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     SW = SWATHobj(DEM,'dx',200,'dy',200);
%     G = gradient8(DEM,'degree');
%     SWG = mapswath(SW,G);
%     plotdz(SWG)
%
% Author: Dirk Scherler (scherler[at]gfz-potsdam.de) and
%         Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: May, 2025

arguments
    SW  SWATHobj
    M
    method = 'linear'
end
   
OUT = SW;
OUT.Z(:) = nan;

if isa(M,'GRIDobj')

    method = validatestring(method,...
        {'linear','nearest','spline','pchip','cubic'},'mapswath','method',3);

    ix = find(~isnan(SW.Z));
    OUT.Z(ix) = interp(M,OUT.X(ix),OUT.Y(ix),method);
    
    % check results
    if any(isnan(OUT.Z))
        warning('Output SWATHobj contains NaN. Input arguments may not fully overlap.')
    elseif all(isnan(OUT.Z))
        warning('Output SWATHobj contains only NaN. Input arguments may not overlap.')
    end

    OUT.name = {'SWATHobj created from SWATHobj:';SW.name};

elseif isa(M,'STREAMobj')
    S = M;

    % Get node attribute list
    if ~isnal(S,method)
        error('Third input must be a node-attribute list')
    else
        nal = method;
    end

    % Find nearest swath pixel for each stream pixel
    [IDX,D] = knnsearch([SW.X(:) SW.Y(:)],[S.x S.y]);
    I = D< max([SW.dx SW.dy])*2;
    OUT.Z(IDX(I)) = nal(I);

    OUT.name = {'SWATHobj created from SWATHobj:';SW.name};

end








