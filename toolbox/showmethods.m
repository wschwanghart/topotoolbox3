function varargout = showmethods(classname,showlink)

%SHOWMETHODS Display class method names and H1 lines in the command line
%
% Syntax
%
%     showmethods(classname)
%     showmethods(classname,showlink)
%
% Description
%
%     showmethods outputs a list of methods available for a specific class
%     (e.g. GRIDobj) in the command window. Only the methods are listed
%     that feature an H1 line, e.g., the first line of the help text block
%     in a function.
%
% Input arguments
%
%     classname     string with class name (e.g. GRIDobj, FLOWobj)
%     showlink      false or true (default). If true, the function displays
%                   hyperlinks to the documentation.
%                   
% Example
%
%     showmethods('GRIDobj')
%
% See also: methods, properties, class
% 
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 2. November, 2024

arguments
    classname 
    showlink (1,1) = true
end

m = methods(classname);

C = cellfun(@(x) numel(x),m,'uniformoutput',true);

if nargout == 1
    IsMakeTable = true;
    OUT = table([],[],'VariableNames',{'Name','Description'});
else
    IsMakeTable = false;
end

maxmethcharacter = max(C);

for r = 1:numel(m)
    s = which([char(classname) '/' m{r}]);
    fileID = fopen(s,'r');
    H1line = false;
    
    try
        while ~H1line && ~feof(fileID)
            tline = fgetl(fileID);
            if isempty(tline)
                continue
            end
            if strcmp(tline(1),'%')
                H1line = true;
            end
        end
        
        if feof(fileID)
            fclose(fileID);
            continue
        end
        
        methodstr = m{r};
        if numel(methodstr)<= maxmethcharacter
            addblanks = repmat(' ',1,maxmethcharacter - numel(methodstr));
        end
        
        h1str = tline(2:end);
        ix    = strfind(h1str,' ');
        h1str = h1str(ix(1)+1:end);

        if ~IsMakeTable
        if ~showlink                
            disp([ upper(methodstr) addblanks ' : ' h1str]);
        else
            disp(['<a href="matlab: doc ' classname '/' methodstr '">' upper(methodstr) '</a>' addblanks ' : ' h1str]);
        end
        else
            h1str(1) = upper(h1str(1));
            t = {string(upper(methodstr)), string(h1str)};
            OUT = [OUT; cell2table(t,'VariableNames',{'Name','Description'})];
        end
                
        fclose(fileID);
        
    catch
        disp(m{r})
    end
    
end
           
if nargout == 1
    varargout{1} = OUT;
end

           


