function Rim = georef2imref2d(varargin)

narginchk(1,2)

if nargin == 1
    % Input is a georeferencing object
    R = varargin{1};
    iscellref = strcmp(R.RasterInterpretation,'cell');

    
else
    % Input is a worldfilematrix and the size of the grid
    wf = varargin{1};
    siz = varargin{2}(1:2);

    % If worldfilematrix, then grid will have cell reference
    iscellref = true;

end
