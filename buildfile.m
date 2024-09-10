function plan = buildfile
    plan = buildplan(localfunctions);
    
    plan("package") = matlab.buildtool.Task( ...
        Description = "Package toolbox", ...
        Dependencies = ["check" "test"], ...
        Actions = @packageToolbox);

    plan.DefaultTasks = ["check" "test" "package"];
end

function checkTask(~)
    issues = codeIssues("toolbox");

    errors = issues.Issues(issues.Issues.Severity=="error",:);
    warnings = issues.Issues(issues.Issues.Severity=="warning",:);

    if ~isempty(errors)
        disp("Errors");
        disp("======")
        disp(errors);
    end

    if ~isempty(warnings)
        disp("Warnings");
        disp("========")
        disp(warnings);
    end

    assert(isempty(errors));
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
