classdef GRIDobj_test < matlab.unittest.TestCase

    methods (Test)
        % Test methods

        function createGRIDobj(test_case)
			% creates a temporary folder and sets it as the current working folder
            test_case.applyFixture(...
                matlab.unittest.fixtures.WorkingFolderFixture);
				
			% Create a GRIDobj from matrices
			for nxy = 3:7
				[X,Y] = meshgrid(1:nxy,1:nxy);
				Z     = rand(size(X));
				DEM   = GRIDobj(X,Y,Z);
				

            inittbx("banana");
            
            directories_to_check = [ ...
                fullfile("banana","toolbox")
                fullfile("banana","toolbox","examples")
                fullfile("banana","tests") ];

            for k = 1:length(directories_to_check)
                d = directories_to_check(k);
                test_case.verifyTrue(directoryExists(d),d);
            end

            files_to_check = [ ...
                fullfile("banana","toolbox","gettingStarted.mlx")
                fullfile("banana","README.md")
                fullfile("banana","LICENSE.md")
                fullfile("banana","buildfile.m")
                fullfile("banana","packageToolbox.m")
                fullfile("banana","toolboxOptions.m") ];

            for k = 1:length(files_to_check)
                f = files_to_check(k);
                test_case.verifyTrue(fileExists(f),f);
            end
        end

        function specifyOutputFolder(test_case)
            test_case.applyFixture(...
                matlab.unittest.fixtures.WorkingFolderFixture); 

            inittbx("banana",OutputFolder = "fruit");

            test_case.verifyTrue(directoryExists(...
                fullfile("fruit","banana","toolbox")));
        end

        function specifyFunctionName(test_case)
            test_case.applyFixture(...
                matlab.unittest.fixtures.WorkingFolderFixture); 

            inittbx("banana",FunctionName = "peel");

            test_case.verifyTrue(fileExists(...
                fullfile("banana","toolbox","peel.m")));
        end

        function specifyToolboxName(test_case)
            test_case.applyFixture(...
                matlab.unittest.fixtures.WorkingFolderFixture); 

            expected = "Fruit Toolbox";
            inittbx("banana",ToolboxName = expected);
            cd("banana")
            opts = toolboxOptions;
            actual = opts.ToolboxName;

            test_case.verifyEqual(actual,expected);
        end

        function specifyToolboxVersion(test_case)
            test_case.applyFixture(...
                matlab.unittest.fixtures.WorkingFolderFixture); 

            expected = "3.1.4";
            inittbx("banana",ToolboxVersion = expected);
            cd("banana")
            opts = toolboxOptions;
            actual = opts.ToolboxVersion;

            test_case.verifyEqual(actual,expected);            
        end

        function runBuild(test_case)
            test_case.applyFixture(...
                matlab.unittest.fixtures.WorkingFolderFixture); 

            inittbx("banana")
            cd banana
            evalc("buildtool");

            test_case.verifyTrue(fileExists(...
                fullfile("release","banana Toolbox.mltbx")));
        end

        function rootNameNotVarname(test_case)
            test_case.applyFixture(...
                matlab.unittest.fixtures.WorkingFolderFixture); 

            inittbx("_banana")

            test_case.verifyTrue(fileExists(...
                fullfile("_banana","toolbox","myfunction.m")));
        end

        function functionNameSpecifiedWithExtension(test_case)
            test_case.applyFixture(...
                matlab.unittest.fixtures.WorkingFolderFixture);  

            inittbx("banana",FunctionName = "peel.m");

            test_case.verifyTrue(fileExists(...
                fullfile("banana","toolbox","peel.m")));
        end
    end
end

function tf = directoryExists(dirname)
    tf = (exist(dirname,"dir") ~= 0);
end

function tf = fileExists(dirname)
    tf = (exist(dirname,"file") ~= 0);
end
