function tt2path
%TT2PATH Add TopoToolbox to MATLAB search path
%
% Syntax
%
%     tt2path
%
% Description
%
%     This function makes it easy to add TopoToolbox to MATLAB's search
%     path.
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 12. October, 2024

% where is tt2path located? 
p = mfilename('fullpath');
f = fileparts(p);
oldfolder = cd(f);

cd("..")
addpath(genpath('toolbox'))

% rmpath(fullfile('toolbox','help'))
% rmpath(fullfile('toolbox','docs','html'))

cd(oldfolder)
