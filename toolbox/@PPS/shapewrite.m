function shapewrite(P,filename,options)

%SHAPEWRITE Write point pattern to shapefile
%
% Syntax
%
%     shapewrite(P,filename)
%     shapewrite(P,filename,pn,pv,...)
%
% Description
%
%     shapewrite writes both the point pattern and the stream network to a
%     shapefile. The points are written to a point shapefile
%     and the stream network to a line shapefile. 
%
% Input arguments
%
%     P         point pattern (PPS)
%     filename  filename of the output shapefile (no *.shp required)
%
%     Parameter name/value pairs
%
%     'streams'  {false} or true. If false, shapewrite will not write the 
%                shapefile of the stream network
%     'postfix'  {'_stream'} postfix for the stream network
%     'marks'    attribute data for each point. This can be either a cell
%                array of size 2xn where n is the number of attributes. The
%                elements in the cell array must contain the name of the
%                attribute and the data. Data can be GRIDobj or node
%                attribute lists. Alternatively, one can supply a table
%                with marks.
%
% Output arguments
%
% See also: STREAMobj/STREAMobj2mapstruct 
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 12. November, 2024

arguments
    P PPS
    filename {mustBeTextScalar}
    options.streams = false
    options.postfix = '_stream'
    options.marks   = []
end

[filepath,name,~] = fileparts(filename);


xy = P.ppxy;
id = (1:npoints(P))';
sh = struct('Geometry','Point',...
            'X',num2cell(xy(:,1)),...
            'Y',num2cell(xy(:,2)),...
            'id',num2cell(id));
if ~isempty(options.marks)
    if iscell(options.marks)
        for r = 1:2:numel(options.marks)
            mark = num2cell(options.marks{r+1});
            [sh.(options.marks{r})] = mark{:};
        end
    elseif istable(options.marks)
        sh = struct2table(sh);
        sh = [sh options.marks];
        sh = table2struct(sh);
%         t = struct2table(options.marks);
%         sh = [sh;t];
    end
end

sh = mapstruct2geotable(sh,"CoordinateReferenceSystem",parseCRS(P),...
    "coordinateSystemType","planar");

shapewrite(sh,fullfile(filepath,name))

if options.streams
    MS = STREAMobj2geotable(P.S);
    shapewrite(MS,fullfile(filepath,name + options.postfix));
end
