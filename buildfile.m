function plan = buildfile
    plan = buildplan(localfunctions);
    
    plan("package") = matlab.buildtool.Task( ...
        Description = "Package toolbox", ...
        Dependencies = ["check" "test"], ...
        Actions = @packageToolbox);

    plan("benchmark").Dependencies = "compile";

    plan.DefaultTasks = ["check" "test" "package"];
end

function checkTask(~)
    issues = codeIssues("toolbox");

    errors = issues.Issues(issues.Issues.Severity=="error", ...
        ["FullFilename", "LineStart", "Description"]);
    warnings = issues.Issues(issues.Issues.Severity=="warning", ...
        ["FullFilename", "LineStart", "Description"]);

    if ~isempty(errors)
        disp("Errors");
        disp("======")
        for i = 1:height(errors)
            fprintf("%s:%d\n",errors{i,"FullFilename"}, ...
                errors{i,"LineStart"});
            fprintf("  %s\n",errors{i,"Description"});
        end
    end

    if ~isempty(warnings)
        disp("Warnings");
        disp("========")
        for i = 1:height(warnings)
            fprintf("%s:%d\n",warnings{i,"FullFilename"}, ...
                warnings{i,"LineStart"});
            fprintf("  %s\n",warnings{i,"Description"});
        end
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

function benchmarkTask(~)
    oldpath = addpath(genpath("toolbox"));
    finalize = onCleanup(@()(path(oldpath)));

    results = runperf("tests/testSnapshot.m");
    disp(sampleSummary(results));
end