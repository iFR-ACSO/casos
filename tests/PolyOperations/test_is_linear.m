classdef test_is_linear < matlab.unittest.TestCase

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

        function test_is_linear_single(testCase,testValue,refValue)
            
            x = casos.Indeterminates('x',2,1);
            % Lyapunov function candidate
            V = casos.PS.sym('v',monomials(x,2:4));
            
            % SOS multiplier
            s2 = casos.PS.sym('s2',monomials(x,0:2:4));
            
            % enforce positivity
            l = 1e-6*(x'*x);
            
          f = [x(1) + x(2);
               x(1)*x(2) - x(1)];
            
           x1 = [V; s2];
        
            % constraints
            g = [s2; 
                 V-l; 
                 s2*(V-1)-nabla(V,x)*f-l];

         for k = 1:3
            actSolution(k) = is_linear(g(k),x1);
         end
            refSolution = logical([1;1;0]');

            
            % Perform assertions if needed
            testCase.verifyEqual(actSolution, refSolution ,"AbsTol",1e-12);

            
        end


    end
   
end
