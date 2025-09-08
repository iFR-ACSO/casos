classdef test_dot < matlab.unittest.TestCase

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
           refValue = [refValuek; refValuek];
        end

         testValue = {testValue};
         refValue  = {refValue};

       end


    end


    methods (Test)

        function test_dot_single(testCase,testValue,refValue)
            
            x = casos.Indeterminates('x',4);
            g = casos.PS.sym('g',monomials(x,2));

            actSolution = full(casadi.DM(full(dot(testValue(1), testValue(1)))));
            
            coeff1 = poly2basis(testValue(1));
            coeff2 = poly2basis(testValue(1));

            refSolution = full(casadi.DM(full(coeff1'*coeff2)));

            % Perform assertions if needed
            testCase.verifyEqual(actSolution, refSolution ,"AbsTol",1e-12);

            
        end


    end
   
end
