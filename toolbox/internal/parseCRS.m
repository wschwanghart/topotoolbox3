function [CRS,isproj] = parseCRS(A)

%PARSECRS Parse coordinate reference system from whatever input
%
% Syntax
%
%     [CRS,isproj] = parseCRS(A)
%
% Description
%
%     The function takes the input A and returns the coordinate reference
%     system and whether the CRS is projected or not. A can be a EPSG ID, a
%     geotable, a GRIDobj or other TopoToolbox object, a projcrs or geocrs.
%
%     If A is a TT object (GRIDobj, FLOWobj, ...), then following 
%     
%     
%     
%
% See also: table2geotable
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 2. September, 2024

if isgeotable(A)
    try
        CRS = A.Shape.GeographicCRS;
        isproj = false;
        return
    catch
        CRS = A.Shape.ProjectedCRS;
        isproj = true;
        return
    end

elseif isTTObject(A)
    
    if isa(A,'PPS')
        A = A.S;
    end

    if isProjected(A)
        CRS = A.georef.ProjectedCRS;
        isproj = true;
        return
    elseif isGeographic(A)
        CRS = A.georef.GeographicCRS;
        isproj = false;
        return
    else
        CRS = [];
        isproj = nan;
        % [~,classstr] = isTTObject(A);
        % error([ classstr ' does not have a defined coordinate reference system.'])
        return
    end

elseif isnumeric(A)
    try 
        CRS = geocrs(A);
        isproj = false;
        return
    catch
        CRS = projcrs(A);
        isproj = true;
        return
    end

elseif isa(A,'geocrs')
    CRS = A;
    isproj = false;
    return

elseif isa(A,'projcrs')
    CRS = A;
    isproj = true;
    return

else
    error('Cannot understand input.')
end
