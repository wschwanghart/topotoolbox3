classdef testHasLibTopoToolbox < matlab.unittest.TestCase
    methods(Test)
        function check_haslibtopotoolbox(testCase)
            verifyReturnsTrue(testCase,@haslibtopotoolbox)
        end
    end
end
