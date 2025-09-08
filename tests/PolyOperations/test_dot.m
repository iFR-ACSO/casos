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
            

            actSolution = full(casadi.DM(full(dot(testValue(1), testValue(1)))));
            
            coeff1 = poly2basis(testValue(1));
            coeff2 = poly2basis(testValue(1));

            refSolution = full(casadi.DM(full(coeff1'*coeff2)));

            % Perform assertions if needed
            testCase.verifyEqual(actSolution, refSolution ,"AbsTol",1e-12);

            
        end

        % function test_dot_multiple(testCase,testValue,refValue)
        % 
        %     actSolution = dot(testValue(1:10), 1);
        % 
        %     refSolution = dot(refValue(1:10), 1);
        % 
        %     c_sopt =  [];
        %     c_cas  =  [];
        %     for k = 1:length(refSolution)
        %             c_sopt = [c_sopt;full(poly2basis(refSolution(k)))];
        %             c_cas  = [c_cas ;full(casadi.DM(poly2basis(actSolution(k))))];
        %     end
        % 
        % 
        %     % Perform assertions if needed
        %     testCase.verifyEqual(c_cas, c_sopt ,"AbsTol",1e-12);
        % end

    end
   
end
