function varargout = gettrace(SW,outformat)
%GETOUTLINE Get trace of swath
%
% Syntax
%
%     out   = gettrace(SW,outformat)
%     [x,y] = gettrace(SW,'xy')
%     [lat,lon] = gettrace(SW,'latlon')
%
% Description
%
%     gettrace returns the trace of a SWATHobj.  
%
% Input arguments
%
%     SW   SWATHobj
%    
%     Parameter name/value pairs
%
%     'output'     'xy','mappolyshape','maplineshape','geopolyshape',
%                  'maplineshape','latlon'
%
% Output arguments
%
%     OUT      Any of the output defined above.
%     x,y      nan-punctuated coordinate vectors
%     lat,lon  nan-punctuated coordinate vectors (latitude, longitude)
%
% Example: See SWATHobj/getoutline
%
% See also: SWATHobj/getoutline
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 12. September, 2024

arguments
    SW  SWATHobj
    outformat = 'maplineshape'
end

if nargout == 1
    validoutformats = {'xy','maplineshape','geopolyshape','latlon'};
else
    validoutformats = {'xy','latlon'};
end
    
outformat = validatestring(outformat,validoutformats,2);

CX{1} = SW.xy(:,1);
CY{1} = SW.xy(:,2);

switch outformat
    case {'geopolyshape','geolineshape','latlon'}
        CRS = parseCRS(SW);
        [CY,CX] = cellfun(@(x,y)projinv(CRS,x,y),CX,CY,'UniformOutput',false);
end

switch outformat
    case {'mappolyshape','maplineshape'}
        outfun = str2func(outformat);
        CX  = cellfun(@(x) x(:)',CX,'UniformOutput',false);
        CY  = cellfun(@(x) x(:)',CY,'UniformOutput',false);
        OUT = outfun(CX,CY);
        OUT.ProjectedCRS = SW.georef.ProjectedCRS;
    case {'geopolyshape','geolineshape'}
        outfun = str2func(outformat);

        CX  = cellfun(@(x) x(:)',CX,'UniformOutput',false);
        CY  = cellfun(@(x) x(:)',CY,'UniformOutput',false);

        OUT = outfun(CY,CX);
        OUT.GeographicCRS = parseCRS(4326);
    case 'xy'
        CX = cellfun(@(x) [x;nan],CX,'UniformOutput',false);
        CY = cellfun(@(x) [x;nan],CY,'UniformOutput',false);
        OUT = [vertcat(CX{:}), vertcat(CY{:})];
    case 'latlon'
        CX = cellfun(@(x) [x;nan],CX,'UniformOutput',false);
        CY = cellfun(@(x) [x;nan],CY,'UniformOutput',false);
        OUT = [vertcat(CY{:}), vertcat(CX{:})];
end

if nargout == 1
    varargout{1} = OUT;
else
    varargout{1} = OUT(:,1);
    varargout{2} = OUT(:,2);
end

