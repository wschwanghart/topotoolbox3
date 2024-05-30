function [cs,ca] = wf2cellsize(wf)

%WF2CELLSIZE Retrieve cellsize from world-file matrix
%
% Syntax
%
%     [cs,ca] = wf2cellsize(wf)
%
% Description
%
%     

% if there is no rotation, wf(2,1) and wf(1,2) will be zero

if wf(2) == 0
    cs = abs(wf(1,1));
    ca = cs^2;
else
    csx = hypot(wf(1,1),wf(2,1));
    csy = hypot(wf(1,2),wf(2,2));
    ca  = csx*csy;
end