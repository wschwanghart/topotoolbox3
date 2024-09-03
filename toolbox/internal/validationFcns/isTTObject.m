function [tf,TTobj] = isTTObject(A)
%isTTObject Determines whether a variable is an instance of a TT object
%
% Syntax
%
%     tf = isTTObject(A)
%     [tf,str] = ...
%
% Description
%
%     isTTObject determines whether a variable is a GRIDobj, FLOWobj,
%     STREAMobj, PPS, DIVIDEobj, or SWATHobj.
%
% Input arguments
%
%     A    any MATLAB variable
%     
% Output arguments
%
%     tf     true, if input is a TT object
%     str    class of the object
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     [tf,str] = isTTObject(DEM)
%     
%     % Should return tf = true and str = 'GRIDobj'
%
% See also: GRIDobj, FLOWobj, STREAMobj
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 2. September, 2024

tf = isa(A,'GRIDobj')  || isa(A,'FLOWobj') || isa(A,'STREAMobj') || ...
     isa(A,'SWATHobj') || isa(A,'PPS')     || isa(A,'DIVIDEobj');

if nargout == 2
    TTobj = class(A);
end