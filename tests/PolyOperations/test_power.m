classdef test_power < matlab.unittest.TestCase

    properties (TestParameter)
        testValue  % test polynmials
        refValue   % reference values
    end
    

    methods (TestParameterDefinition,Static)

       function [testValue,refValue] = initializeTestData()

           testValue = [];
           refValue  = [];

        for k = 1:5
           [testValuek,refValuek] = genTestPoly();
           testValue = [testValue; testValuek];
           refValue = [refValue; refValuek];
        end

         testValue = {testValue};
         refValue  = {refValue};

       end


    end


    methods (Test)

        function test_times_power(testCase,testValue,refValue)

            actSolution = power(testValue(1), 3);

            refSolution = power(refValue(1), 3);

            c_sopt = full(poly2basis(refSolution));
            c_cas  = full(casadi.DM(poly2basis(actSolution)));
            
            % Perform assertions if needed
            testCase.verifyEqual(c_cas, c_sopt ,"AbsTol",1e-12);

            
        end

        % function test_power_multiple(testCase,testValue,refValue)
        % 
        %     actSolution = power(testValue(1:2),5);
        % 
        %     refSolution = power(refValue(1:2), 5);
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
