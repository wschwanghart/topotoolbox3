function [V,IX] = mapfromnal(FD,S,nal,cl)

%MAPFROMNAL Map values from node-attribute list to nearest upstream grid
%
% Syntax
%
%     V = mapfromnal(FD,S,nal)
%     [V,IX] = mapfromnal(FD,S,nal)
%
% Description
%
%     mapfromnal takes a STREAMobj S and an associated node-attribute list
%     nal and maps the values in the nal to the nearest grid values
%     measured along flowpaths based on the FLOWobj FD. S should be have
%     been derived from FD.
%
%     If FD is a multiple flow direction object and S and FD have been
%     created with the function multi2single, then stream values cannot be
%     unambiguously mapped to upstream pixels. In this case, values
%     represent averages of the values of stream pixels which hillslopes 
%     are contributing to. 
%
%     
%
% Input arguments
%
%     FD      FLOWobj
%     S       STREAMobj
%     nal     node-attribute list
%
% Output arguments
%
%     V       GRIDobj with values derived from nal
%     IX      matrix with size V.size with linear indices into 
%             node-attributes of S. Elements in IX with no downstream
%             stream pixel are zero. IX is empty, if FD is a multiple flow
%             direction object.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     A = flowacc(FD);
%     V = mapfromnal(FD,S,A);
%     imagesc(V)
%
% See also: FLOWobj, STREAMobj, flowdistance, vertdistance2stream
% 
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 31. August, 2024

arguments
    FD    FLOWobj
    S     STREAMobj
    nal   {mustBeGRIDobjOrNal(nal,S)}
    cl = 'single'
end


switch lower(FD.type)
    case 'single'
        % If the variable nal is a GRIDobj than extract the nal
        nal = ezgetnal(S,nal);

        % Use propagatevaluesupstream
        IX = propagatevaluesupstream(FD,S.IXgrid,uint32(1:numel(S.IXgrid)),...
            'fillval',zeros(1,'uint32'),'overwrite',false);
        IX = IX.Z;
        V  = GRIDobj(FD,nan(FD.size,cl));
        I  = IX>0;
        V.Z(I) = nal(IX(I));

    otherwise

        DEMS = STREAMobj2GRIDobj(S,nal);
        DEMS.Z(isnan(DEMS.Z)) = 0;
        ISSTREAM = STREAMobj2GRIDobj(S);

        % There are some pixels that do not have downstream stream pixels. Usually,
        % these are part of smaller catchment along the DEM edges. We remove links
        % from the flow network that connect any of these pixels.
        I   = dependencemap(FD,S.IXgrid);
        ii   = I.Z(FD.ix) & I.Z(FD.ixc);
        FD.ix = FD.ix(ii);
        FD.ixc = FD.ixc(ii);

        % Accumulated fractions
        FR  = zeros(DEMS.size);

        for r = numel(FD.ix):-1:1
            FR(FD.ix(r)) = FD.fraction(r) + FR(FD.ix(r));
        end

        for r = numel(FD.ix):-1:1
            if ISSTREAM.Z(FD.ix(r))
                continue
            end

            DEMS.Z(FD.ix(r)) = FD.fraction(r).*DEMS.Z(FD.ixc(r)) .* 1/FR(FD.ix(r)) ...
                + DEMS.Z(FD.ix(r));
        end

        V = DEMS;
        % V = clip(V,~ISSTREAM);
        V.Z(~I.Z) = nan;
        if nargout == 2
            IX = [];
            warning('Second output is empty if FD contains multiple flow directions.')
        end
end
