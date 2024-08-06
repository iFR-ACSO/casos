classdef test_nabla < matlab.unittest.TestCase

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

        function test_nabla_single(testCase,testValue,refValue)
            
            actSolution = nabla(testValue(1), testValue(1).indeterminates);
      
            refSolution = jacobian(refValue(1), mpvar('x',refValue(1).nvar,1));

            c_sopt = full(poly2basis(refSolution));
            c_cas  = full(casadi.DM(poly2basis(actSolution)));
            
            % Perform assertions if needed
            testCase.verifyEqual(c_cas, c_sopt(:) ,"AbsTol",1e-12);

            
        end

   

    end
   
end
