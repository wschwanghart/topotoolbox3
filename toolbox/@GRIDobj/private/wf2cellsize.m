function cs = wf2cellsize(wf)

% This function extracts the cellsize from a world-file matrix

if wf(2,1) == 0
    if abs(wf(1,1)) ~= abs(wf(2,2))
        warning("Cells are not square.")
    end
    cs = abs(wf(1));
else
    cs1 = hypot(wf(1,1),wf(2,1));
    cs2 = hypot(wf(1,2),wf(2,2));

    if cs1 == cs2
        cs = cs1;
    elseif abs(cs1-cs2) < eps*1000
        cs = cs1;
    else
        cs = cs1;
        warning("Cells are not square.")
    end
end
