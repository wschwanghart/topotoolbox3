function cmap = loadcmap(cmaptype)

%LOADCMAP load colormap

p = fileparts(mfilename('fullpath'));
p = char(p);
cmap = load([p filesep 'colormaps' filesep cmaptype '.mat']);
cmap = cmap.(cmaptype);