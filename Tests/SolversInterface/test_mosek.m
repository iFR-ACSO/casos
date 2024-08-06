classdef test_mosek < matlab.unittest.TestCase
    % WARNING: THIS IS JUST  A PLACEHOLDER FILE

    properties (TestParameter)
        testValue  % test polynmials
        refValue   % reference values
    end
    

    methods (TestParameterDefinition,Static)

       function [testValue,refValue] = initializeTestData()

           testValue = {1}; 
           refValue  = {1};

       end
    end


    methods (Test)

        function test_dot_single(testCase,testValue,refValue)
            

            actSolution = 1; %full(casadi.DM(full(dot(testValue(1), testValue(1)))));
            
            % coeff1 = poly2basis(testValue(1));
            % coeff2 = poly2basis(testValue(1));

            refSolution = 2; %full(casadi.DM(full(coeff1'*coeff2)));

            % Perform assertions if needed
            testCase.verifyEqual(actSolution, refSolution ,"AbsTol",1e-12);

            
        end

    end
   
end