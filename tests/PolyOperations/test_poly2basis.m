classdef test_poly2basis < matlab.unittest.TestCase

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

        function test_poly2basis_single(testCase,testValue,refValue)

            c_sopt = full(poly2basis(refValue(1)));
            c_cas  = full(casadi.DM(poly2basis(testValue(1))));
            
            % Perform assertions if needed
            testCase.verifyEqual(c_cas, c_sopt ,"AbsTol",1e-12);

            
        end

        function test_poly2basis_multiple(testCase,testValue,refValue)

            
            c_sopt =  [];
            c_cas  =  [];
            for k = 1:length(refValue)
                    
                    c_sopt = [c_sopt;full(poly2basis(refValue(k)))];
                    c_cas  = [c_cas ;full(casadi.DM(poly2basis(testValue(k))))];
            end
  

            % Perform assertions if needed
            testCase.verifyEqual(c_cas, c_sopt ,"AbsTol",1e-12);
        end

    end
   
end
