classdef testSnapshot < matlab.unittest.TestCase
    % testSnapshot Snapshot testing using topotoolbox3 as a reference
    % The data for the snapshot tests are stored in the
    % TopoToolbox/test_data repository. This repository is added as a
    % submodule to topotoolbox3 at the path tests/data.
    %
    % The basic strategy for these tests:
    %
    % 1. The snapshot tests are run just like a normal test for
    %    topotoolbox3.
    % 2. If there is no data for a given test in the tests/data directory,
    %    the test will create it.
    % 3. If the appropriate data does exist in the tests/data directory,
    %    they will be compared against the results computed here. The test
    %    will fail if the saved version and the computed version do not
    %    match according to the test. Be careful when adding new tests
    %    requiring floating point comparison.
    % 4. Snapshotting of results must currently be done manually. If you
    %    want to save the results of a test run:
    %    - Commit the new data to the /submodule/, not to topoboolbox3
    %    - Push the commit to your fork of TopoToolbox/test_data and make a
    %    pull request.
    %    - Once the pull request is merged, pull the new changes into the
    %    submodule and commit the changed submodule to the
    %    topotoolbox3 repository.
    % 5. If a change to the results is known and expected, delete the
    %    appropriate snapshot from the tests/data repository and create the
    %    new one following the procedure outlined above.

    properties
        dem
    end

    properties (ClassSetupParameter)
        dataset
    end

    methods (TestParameterDefinition,Static)
        function dataset = findDatasets()
            % Find all the existing snapshot datasets
            available_datasets = [{},struct2table(dir("snapshots/*/dem.tif")).folder];
            if ~isempty(available_datasets)
                dataset = available_datasets;
            else
                error("No snapshots found.");
            end
        end
    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function read_data(testCase,dataset)
            data_file = fullfile(dataset,"dem.tif");
            testCase.dem = GRIDobj(data_file);
        end
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        % Test methods
        function fillsinks(testCase,dataset)
            demf = testCase.dem.fillsinks();

            result_file = fullfile(dataset,"fillsinks.tif");

            if ~isfile(result_file)
                % Write the result to the directory
                demf.GRIDobj2geotiff(result_file);
            else
                demf_result = GRIDobj(result_file);
                testCase.verifyEqual(demf_result.Z,demf.Z);
            end
        end

        function identifyflats(testCase, dataset)
            demf = testCase.dem.fillsinks();
            [FLATS, SILLS, CLOSED] = demf.identifyflats();

            flats_file = fullfile(dataset,"identifyflats_flats.tif");
            sills_file = fullfile(dataset,"identifyflats_sills.tif");
            closed_file = fullfile(dataset,"identifyflats_closed.tif");

            if ~isfile(flats_file)
                FLATS.GRIDobj2geotiff(flats_file);
            else
                flats_result = GRIDobj(flats_file);
                testCase.verifyEqual(logical(flats_result.Z),FLATS.Z);
            end

            if ~isfile(sills_file)
                SILLS.GRIDobj2geotiff(sills_file);
            else
                sills_result = GRIDobj(sills_file);
                testCase.verifyEqual(logical(sills_result.Z),SILLS.Z);
            end

            if ~isfile(closed_file)
                CLOSED.GRIDobj2geotiff(closed_file);
            else
                closed_result = GRIDobj(closed_file);
                testCase.verifyEqual(logical(closed_result.Z),CLOSED.Z);
            end
        end
    end
end