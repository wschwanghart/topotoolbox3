function plan = buildfile
    plan = buildplan(localfunctions);
    
    plan("check") = matlab.buildtool.tasks.CodeIssuesTask("toolbox");
    plan("package") = matlab.buildtool.Task( ...
        Description = "Package toolbox", ...
        Dependencies = ["check" "test"], ...
        Actions = @packageToolbox);

    plan.DefaultTasks = ["check" "test" "package"];
end

function testTask(~)
    oldpath = addpath(genpath("toolbox"));
    finalize = onCleanup(@()(path(oldpath)));

    results = runtests("IncludeSubfolders",true);
    disp(results);
    assertSuccess(results);
end

function compileTask(~)
    %% Path to MATLAB provided cmake
    cmake = fullfile(matlabroot, 'bin', computer('arch'),'cmake','bin','cmake');
    % The windows version has a blank in the program folder name
    % (C:/Program Files), thus we use additional double quotes.
    cmake = ['"' cmake '"'];

    if system(strcat(cmake," --version")) ~= 0
        disp("CMake is not provided by MATLAB. MATLAB Coder must be installed.");
        return;
    end
    %% Run cmake
    system(strcat(cmake," -S bindings/ -B build -DCMAKE_BUILD_TYPE=Release"));
    system(strcat(cmake," --build build --config Release"));
    system(strcat(cmake," --install build --prefix toolbox/internal/mex --component Runtime"));
end
