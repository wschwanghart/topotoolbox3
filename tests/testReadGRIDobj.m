classdef testReadGRIDobj < matlab.unittest.TestCase

    properties      
        url = 'https://raw.githubusercontent.com/TopoToolbox/DEMs/master/';
        currentFolder = pwd
        tempFolder = fullfile(tempdir,'tt3_testdir')
    end
    properties (TestParameter)
        files = {'kedarnath','kunashiri','perfectworld','taalvolcano','tibet'};
    end

    methods (TestClassSetup)
        function readFiles(testCase) 

            % Create a folder in the temporary directory
            if ~exist(testCase.tempFolder,"file")
                mkdir(testCase.tempFolder)
            end
            
            % add path with tests to path
            p = mfilename('fullpath');
            [p,~] = fileparts(p);
            addpath(genpath(p))

            testCase.currentFolder = cd(testCase.tempFolder);
            % Now we should be in the test folder

            
            % Check whether the test files are already downloaded
            for r = 1:numel(testCase.files)
                % provide absolute paths
                fn = [testCase.files{r} '.tif'];
                if ~(exist(fullfile(testCase.tempFolder,fn),'file') == 2)
                    websave(fullfile(testCase.tempFolder,fn),fullfile(testCase.url,fn));
                end
            end 
        end
    end


       
    methods (Test, ParameterCombination = 'sequential')
        % Test methods

        function read_GRIDobj(testCase,files)
				
            % Read example DEM
            fn  = fullfile(testCase.tempFolder,[files '.tif']);
            DEM = GRIDobj(fn);

            verifyInstanceOf(testCase,DEM,'GRIDobj')
        end

        function write_GRIDobj_test1(testCase,files)
            
            % Read GRIDobj
            fn  = fullfile(testCase.tempFolder,[files '.tif']);
            DEM = GRIDobj(fn);
            % Write DEM
            GRIDobj2geotiff(DEM,'test.tif')
            I = exist('test.tif','file');
            verifyEqual(testCase,I,2)

        end

        function write_GRIDobj_test2(testCase,files)

            % Read example DEM
            fn  = fullfile(testCase.tempFolder,[files '.tif']);
            DEM = GRIDobj(fn);
            % Write DEM
            GRIDobj2geotiff(DEM,'test.tif')
            DEM2 = GRIDobj("test.tif");            
            verifyEqual(testCase,DEM.Z,DEM2.Z)
        end

        function write_GRIDobj_test3(testCase,files)

            % Read example DEM
            fn  = fullfile(testCase.tempFolder,[files '.tif']);
            DEM = GRIDobj(fn);

            % Write DEM
            GRIDobj2geotiff(DEM,'test.tif')
            DEM2 = GRIDobj("test.tif");            
            verifyEqual(testCase,DEM.wf,DEM2.wf)
        end

        function write_GRIDobj_test4(testCase,files)

            % Read example DEM
            fn  = fullfile(testCase.tempFolder,[files '.tif']);
            DEM = GRIDobj(fn);
            % Write DEM
            GRIDobj2geotiff(DEM,'test.tif')
            DEM2 = GRIDobj("test.tif");     
            if isProjected(DEM)
                verifyEqual(testCase,DEM.georef.ProjectedCRS,DEM2.georef.ProjectedCRS)
            elseif isGeographic(DEM)
                verifyEqual(testCase,DEM.georef.GeographicCRS,DEM2.georef.GeographicCRS)
            else
                verifyEqual(testCase,DEM.georef,DEM2.georef)
            end

        end

        function projectGRIDobj(testCase)
				
            % Read example DEM
            DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
            DEM2 = project(DEM,4326,'res',0.01);

            verifyEqual(testCase,abs(DEM2.wf(1)),0.01,'AbsTol',1e-6)
            
        end
    end

end
