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
%     geotable, a GRIDobj, a projcrs or geocrs.
%
% See also: table2geotable
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 9. July, 2024

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

elseif isa(A,'GRIDobj') 
    if isProjected(A)
        CRS = A.georef.ProjectedCRS;
        isproj = true;
        return
    elseif isGeographic(A)
        CRS = A.georef.GeographicCRS;
        isproj = false;
        return
    else
        error('GRIDobj does not have a defined coordinate reference system.')
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
