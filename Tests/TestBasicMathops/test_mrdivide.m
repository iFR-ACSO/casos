classdef test_mrdivide < matlab.unittest.TestCase

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

        function test_mrdivide_single(testCase,testValue,refValue)
            
            B = 2;
            actSolution = mrdivide(testValue(1), B);

            refSolution = mrdivide(refValue(1), B);

            c_sopt = full(poly2basis(refSolution));
            c_cas  = full(casadi.DM(poly2basis(actSolution)));
            
            % Perform assertions if needed
            testCase.verifyEqual(c_cas, c_sopt ,"AbsTol",1e-12);

            
        end

        function test_mrdivide_multiple(testCase,testValue,refValue)
        
            actSolution = mrdivide(testValue(1:10), 3);

            refSolution = mrdivide(refValue(1:10), 3);
            
            c_sopt =  [];
            c_cas  =  [];
            for k = 1:length(refSolution)
                    c_sopt = [c_sopt;full(poly2basis(refSolution(k)))];
                    c_cas  = [c_cas ;full(casadi.DM(poly2basis(actSolution(k))))];
            end
  

            % Perform assertions if needed
            testCase.verifyEqual(c_cas, c_sopt ,"AbsTol",1e-12);
        end

    end
   
end
