function [OUT,OUT2] = streampoi(S,type,outformat)

%STREAMPOI Points of interest on the stream network
%
% Syntax
%
%     V = streampoi(S,type)
%     V = streampoi(S,type,outformat)
%     [x,y] = streampoi(S,type,'xy')
%
% Description
%
%     streampoi returns points of interest in a stream network. These
%     points can be 'channelheads', 'confluences' or 'outlets'.
%
% Input arguments
%   
%     S           stream network (class: STREAMobj)
%     type        'channelheads' (default), 'confluences', 'bconfluences' 
%                 or 'outlets'. bconfluences returns the stream pixels that 
%                 are located immediately upstream to confluences.
%                 type can also be a cell array of strings with the
%                 different types listed above.
%     outformat   'xy': nx2 coordinate matrix with x and y coordinates of 
%                       n points (default). With two output arguments, the
%                       function will return two vectors with x and y
%                       coordinates, respectively.
%                 'ix': nx1 vector with linear indices into an instance of
%                       GRIDobj with the same dimension as the GRIDobj from
%                       which S was derived.
%                 'logical': node attribute list (logical) 
%                 'nal': same as 'logical'
%                 'mappoint': mappoint (see function mapshape)
%                 'geopoint': geopoint (see function geopoint). Requires a
%                       valid coordinate reference system.
%                 'PPS' instance of PPS
%
% Output
%
%     V           output as specified in outformat
%     x,y         coordinate vectors
%
% Example 1
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD = FLOWobj(DEM);
%     S  = STREAMobj(FD,'minarea',1000);
%     xy = streampoi(S,'confluence');
%     plot(S);
%     hold on
%     plot(xy(:,1),xy(:,2),'s')
%     hold off
%
% Example 2
%
%     P  = streampoi(S,'confluence','PPS');
%     P.z = DEM;
%     plotdz(P)
% 
% See also: STREAMobj, FLOWobj/streampoi
% 
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 8. November, 2024

arguments
    S   STREAMobj
    type = 'channelheads'
    outformat  {mustBeMember(outformat,...
        {'xy','ix','logical','nal','mappoint','geopoint','PPS'})} = 'xy'
end

% Place type into a cell array 
if ischar(type) || isstring(type)
    type = {type};
end

% Sparse matrix representation of the network
ix  = S.ix;
ixc = S.ixc;
nrc = numel(S.x);
M   = sparse(double(ix),double(ixc),true,nrc,nrc);

% Output logical array
V = false(nrc,1);

% Go through types of points of interest
for r = 1:numel(type)
    t = validatestring(type{r},...
        {'channelheads','confluences','bconfluences','outlets'},...
        'streampoi','type',2);
   
    switch t
        case 'channelheads'
            
            V = (sum(M,1)'==0) & (sum(M,2) ~= 0) | V;
            
        case 'confluences'
            
            V = sum(M,1)'>1 | V;
            
        case 'bconfluences'
            
            V2 = sum(M,1)'>1;
            V2 = any(M(:,V2),2);
            V  = V | V2;
            
        case 'outlets'
            
            V = ((sum(M,2)==0) & (sum(M,1)'~=0)) | V;
    end
end

% prepare output
switch lower(outformat)
    case 'xy'      
        ix   = find(V);
        ix   = ix(:);
        OUT  = [S.x(ix) S.y(ix)];
        if nargout == 2
            OUT2 = OUT(:,2);
            OUT = OUT(:,1);
        end
        
    case 'ix'
        ix   = find(V);
        ix   = ix(:);
        OUT  = S.IXgrid(ix);
        
    case {'logical','nal'}
        
        OUT  = full(V(:));
        
    case 'mappoint'
        
        ix   = find(V);
        ix   = ix(:);
        OUT  = mappoint(S.x(ix),S.y(ix));
        
    case 'geopoint'
        
        ix   = find(V);
        ix   = ix(:);
        [lat,lon] = projinv(S.georef.ProjectedCRS,S.x(ix),S.y(ix));
        OUT  = geopoint(lat,lon);
        
    case 'pps'
        
        OUT   = PPS(S,'pp',S.IXgrid(V));
end