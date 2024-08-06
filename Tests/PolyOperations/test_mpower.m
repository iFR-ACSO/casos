classdef test_mpower < matlab.unittest.TestCase

    properties (TestParameter)
        testValue  % test polynmials
        refValue   % reference values
    end
    

    methods (TestParameterDefinition,Static)

       function [testValue,refValue] = initializeTestData()

           testValue = [];
           refValue  = [];

        for k = 1:10
           [testValuek,refValuek] = genTestPoly();
           testValue = [testValue; testValuek];
           refValue = [refValue; refValuek];
        end

         testValue = {testValue};
         refValue  = {refValue};

       end


    end


    methods (Test)

        function test_mpower_single(testCase,testValue,refValue)

            actSolution = mpower(testValue(1), 2);

            refSolution = mpower(refValue(1), 2);

            c_sopt = full(poly2basis(refSolution));
            c_cas  = full(casadi.DM(poly2basis(actSolution)));
            
            % Perform assertions if needed
            testCase.verifyEqual(c_cas, c_sopt ,"AbsTol",1e-12);

            
        end



    end
   
end
