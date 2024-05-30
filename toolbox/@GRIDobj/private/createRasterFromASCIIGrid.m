function [Z,R,wf] = createRasterFromASCIIGrid(filename,options)

%CREATERASTERFROMASCIIGRID Creates a grid from an Esri ASCII grid
%
% Syntax
%
%     [Z,R] = createRasterFromASCIIGrid(filename,options)
%
% 

arguments
    filename 
    options.OutputType = 'single'
end

% Read raster using readgeoraster (if mapping toolbox available)
if license('test','MAP_Toolbox')
    [Z,R] = readgeoraster(filename,...
        'OutputType',options.OutputType,...
        'Bands',1,...
        'CoordinateSystemType',options.CoordinateSystemType);
    wf = worldFileMatrix(R);
    
    % check, if R contains a projcrs or a geocrs, if any

else
    % if mapping toolbox is not available, use rasterread
    [Z,wf] = rasterread(filename);    
    Z  = cast(Z,options.OutputType);
    R  = [];
end

