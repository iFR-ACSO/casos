classdef test_sedumi < matlab.unittest.TestCase
    % WARNING: THIS IS JUST  A PLACEHOLDER FILE
    
    properties (TestParameter)
        testValue  % test polynmials
        refValue   % reference values
    end
    

    methods (TestParameterDefinition,Static)

       function [testValue,refValue] = initializeTestData()
            
            % verify if sedumi is in path
            if ~exist('sedumi','file')
                error('sedumi is not in path');
            end
            
            % verify if sosopt and multipoly is in path
            if ~exist('sosopt','file') || ~exist('multipoly','file')
                error('sosopt or multipoly is not in path');
            end
            
            % if yes then run everything nicely, otherwise output a flag
            % not to run but should tell that the test didn't run due to
            % missing packages in path

           testValue = {1}; 
           refValue  = {1};

       end

    end


    methods (Test)

        function test_sedumi_sdp(testCase,testValue,refValue)
            
            actSolution = 1; % full(casadi.DM(full(dot(testValue(1), testValue(1)))));
            
            % coeff1 = poly2basis(testValue(1));
            % coeff2 = poly2basis(testValue(1));

            refSolution = 1; %full(casadi.DM(full(coeff1'*coeff2)));

            % Perform assertions if needed
            testCase.verifyEqual(actSolution, refSolution ,"AbsTol",1e-12);

        end

    end
   
end
