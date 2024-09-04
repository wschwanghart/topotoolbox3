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
