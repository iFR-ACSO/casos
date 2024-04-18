classdef test_plus < matlab.unittest.TestCase

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
           refValue  = [refValue; refValuek];
        end

         testValue = {testValue};
         refValue  = {refValue};

       end


    end




    methods (Test)
                %% ------------------------------------------------------------
                %
                %  Input permutation of nominal tests
                %  0 : constant
                %  1 : single polynomial
                %  n : n polynomials
                %         
                %         A     B     C     D     E     F     G     H   I
                %
                %   a:    0     1     n     0     1     n     0     1   n
                %   b:    0     0     0     1     1     1     n     n   n
                %
                % -------------------------------------------------------------
            %% case A
            function testPlus_case_A(testCase,testValue,refValue)
    
                actSolution = plus(casos.PS(5), casos.PS(5));
    
                refSolution = 10;
                c_cas  = full(casadi.DM(poly2basis(actSolution)));
                
                % Perform assertions if needed
                testCase.verifyEqual(c_cas, refSolution ,"AbsTol",1e-12);
    
       
            end
    
            %% case B
            function testPlus_case_B(testCase,testValue,refValue)
    
                actSolution = plus(testValue(1), casos.PS(3));
    
                refSolution = plus(refValue(1), 3);
    
                c_sopt = full(poly2basis(refSolution));
                c_cas  = full(casadi.DM(poly2basis(actSolution)));
                
                % Perform assertions if needed
                testCase.verifyEqual(c_cas, c_sopt ,"AbsTol",1e-12);
    
       
            end
    
            %% case C
            function testPlus_case_C(testCase,testValue,refValue)
    
                actSolution = plus(testValue(1:10), 4);
    
                refSolution = plus(refValue(1:10), 4);
    
                c_sopt =  [];
                c_cas  =  [];
                for k = 1:length(refSolution)
                        c_sopt = [c_sopt;full(poly2basis(refSolution(k)))];
                        c_cas  = [c_cas ;full(casadi.DM(poly2basis(actSolution(k))))];
                end
      
                
                % Perform assertions if needed
                testCase.verifyEqual(c_cas, c_sopt ,"AbsTol",1e-12);
    
       
            end
    
            %% case D
            function testPlus_case_D(testCase,testValue,refValue)
    
                actSolution = plus(4, testValue(1));
    
                refSolution = plus(4,refValue(1));
    
                c_sopt =  [];
                c_cas  =  [];
                for k = 1:length(refSolution)
                        c_sopt = [c_sopt;full(poly2basis(refSolution(k)))];
                        c_cas  = [c_cas ;full(casadi.DM(poly2basis(actSolution(k))))];
                end
      
                
                % Perform assertions if needed
                testCase.verifyEqual(c_cas, c_sopt ,"AbsTol",1e-12);
    
       
            end
    
            %% case E
            function testPlus_case_E(testCase,testValue,refValue)
    
                actSolution = plus(testValue(1), testValue(2));
    
                refSolution = plus(refValue(1), refValue(2));
    
                c_sopt = full(poly2basis(refSolution));
                c_cas  = full(casadi.DM(poly2basis(actSolution)));
                
                % Perform assertions if needed
                testCase.verifyEqual(c_cas, c_sopt ,"AbsTol",1e-12);
    
       
            end
    
    
            %% case F
            function testPlus_case_F(testCase,testValue,refValue)
    
                actSolution = plus(testValue(1:10),testValue(1));
    
                refSolution = plus(refValue(1:10),refValue(1));
    
                c_sopt =  [];
                c_cas  =  [];
                for k = 1:length(refSolution)
                        c_sopt = [c_sopt;full(poly2basis(refSolution(k)))];
                        c_cas  = [c_cas ;full(casadi.DM(poly2basis(actSolution(k))))];
                end
      
                
                % Perform assertions if needed
                testCase.verifyEqual(c_cas, c_sopt ,"AbsTol",1e-12);
    
       
            end
    
            %% case G
            function testPlus_case_G(testCase,testValue,refValue)
    
                actSolution = plus(5,testValue(1:10));
    
                refSolution = plus(5,refValue(1:10));
    
                c_sopt =  [];
                c_cas  =  [];
                for k = 1:length(refSolution)
                        c_sopt = [c_sopt;full(poly2basis(refSolution(k)))];
                        c_cas  = [c_cas ;full(casadi.DM(poly2basis(actSolution(k))))];
                end
      
                
                % Perform assertions if needed
                testCase.verifyEqual(c_cas, c_sopt ,"AbsTol",1e-12);
    
       
            end
    
    
            %% case H
            function testPlus_case_H(testCase,testValue,refValue)
    
                actSolution = plus(testValue(1),testValue(1:10));
    
                refSolution = plus(refValue(1),refValue(1:10));
    
                c_sopt =  [];
                c_cas  =  [];
                for k = 1:length(refSolution)
                        c_sopt = [c_sopt;full(poly2basis(refSolution(k)))];
                        c_cas  = [c_cas ;full(casadi.DM(poly2basis(actSolution(k))))];
                end
      
                
                % Perform assertions if needed
                testCase.verifyEqual(c_cas, c_sopt ,"AbsTol",1e-12);
    
       
            end
    
             %% case I
             function testPlus_case_I(testCase,testValue,refValue)
    
                actSolution = plus(testValue(1:10),testValue(11:20));
    
                refSolution = plus(refValue(1:10),refValue(11:20));
    
                c_sopt =  [];
                c_cas  =  [];
                for k = 1:length(refSolution)
                        c_sopt = [c_sopt;full(poly2basis(refSolution(k)))];
                        c_cas  = [c_cas ;full(casadi.DM(poly2basis(actSolution(k))))];
                end
      
                
                % Perform assertions if needed
                testCase.verifyEqual(c_cas, c_sopt ,"AbsTol",1e-12);
    
       
            end
            
           %% ------------------------------------------------------------
           %
           % Test error cases
           %
           % ------------------------------------------------------------- 
           % function testPlus_case_error(testCase,testValue,refValue)
           % 
           % 
           %      sza = size(testValue(1:10));
           %      szb = size([]);
           %      str1 = sprintf('%dx%d', sza(1), sza(2));
           %      str2 = sprintf('%dx%d', szb(1), szb(2));
           % 
           % 
           % 
           %      errortxt = ['Polynomials have incompatible sizes for this operation','([' str1 '] vs. [' str2 ']).'];
           % 
           %      % Perform assertions if needed
           %      % testCase.verifyEqual(c_cas, c_sopt ,"AbsTol",1e-12);
           %      testCase.verifyError(plus(testValue(1:10),[]),errortxt)
           % 
           %  end

        end

  
end
