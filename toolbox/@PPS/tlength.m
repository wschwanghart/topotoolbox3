function l = tlength(P)

%TLENGTH Total length of the stream network
%
% Syntax
%
%     l = tlength(P)
%
% Description
%
%     tlength returns the total length of the stream network measured in
%     units of distance.
%
% 
% See also: PPS, PPS/npoints 
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 5. July, 2024

arguments 
    P    PPS
end

l = info(P.S,'totallength');