function out = haslibtopotoolbox
%HASLIBTOPOTOOLBOX Test if the libtopotoolbox MEX bindings are available
mexx = exist(['mex/tt_has_topotoolbox.' mexext],'file');
if mexx == 3
    out = tt_has_topotoolbox;
else
    out = false;
end
end
