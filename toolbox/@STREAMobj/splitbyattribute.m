function [CS,zagg,edges] = splitbyattribute(S,z,n)
%SPLITBYATTRIBUTE Create a cell array of STREAMobjs using an attribute
%
% Syntax
%
%     [CS,zagg] = splitbyattribute(S,z,n)
%     [CS,zagg] = splitbyattribute(S,z,edges)
%     [CS,zagg,edges] = ...
%
% Description
%
%     The function takes a STREAMobj S and a numeric node attribute list
%     (or GRIDobj) z and distributes S into a cell array of STREAMobjs so 
%     that the each element contains a STREAMobj which covers one of n
%     partions of z. 
%
%     If you have a categorical node-attribute list z, then use
%     STREAMobj2cell(S,'label',z).
%
% Input arguments
%
%     S      STREAMobj
%     z      numeric node-attribute list or GRIDobj
%     n      number of partitions (default = 10)
%     edges  vector that determines the bin edges (see histcounts)
% 
% Output arguments
%
%     CS     cell array of STREAMobj
%     zagg   central values of the partitioned values of z
%     edges  bin edges (zagg is calculated by
%            (edges(1:end-1)+edges(2:end))/2
%     
% Example 1
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     S = STREAMobj(FD,'minarea',1000);
%     [CS,zagg] = splitbyattribute(S,DEM,10);
%     hold on
%     cellfun(@plot,CS);
%     hold off
%
% See also: STREAMobj/STREAMobj2cell, STREAMobj/wmplot, histcounts
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 18. July, 2024

arguments
    S  STREAMobj
    z  {mustBeGRIDobjOrNal(z,S)}
    n  = 10
end

z = ezgetnal(S,z);

[N,edges,bin] = histcounts(z,n);
zagg = (edges(1:end-1) + edges(2:end))/2;
zagg = zagg(:);

CS = STREAMobj2cell(S,'labels',bin);
zagg = zagg(N>0);




