% ========================================================================
%
% Test Name: test_sdpsol.m
%
% Test Description: 
%   > Check for a simple SDP with only linear and psd cones if the
%   interface with mosek and sedumi works.
%   > see https://yalmip.github.io/tutorial/semidefiniteprogramming/
%               
% Test Procedure: 
%
%   min       c'x       |   c   = [0; ...; 0]       
%   s.t.    x in K_x    |   K_x = single psd cone
%          Ax in K_c    |   K_c = single psd cone 
%                               
% Date: 09/06/2024 
%
% ========================================================================

classdef test_sdpsol_sedumi < matlab.unittest.TestCase

    properties (TestParameter)
        sdp  
        opts 
    end

    properties (Access = private, Constant)
        PackagesAvailable   = test_sdpsol_sedumi.checkRequiredPackages(1);
        MissingPackages     = test_sdpsol_sedumi.checkRequiredPackages(2);
    end
    
    methods (Static)
        function output = checkRequiredPackages(out_select)
            % Define a list of required packages and their check functions
            packages = {'sedumi'};
            packageChecks = @(pkgName) exist(pkgName, 'file') == 2;

            % Initialize the flag for package availability
            available = true;
            missingPackages = {};

            % Check each package
            for i = 1:numel(packages)
                pkgName = packages{i};
                if ~packageChecks(pkgName)
                    available = false;
                    missingPackages{end+1} = pkgName;
                end
            end

            % set output
            if out_select==1
                output = available;
            else
                output = missingPackages;
            end
        end
    end
    
    methods (TestClassSetup)
        function setupClass(testCase)
            if ~test_sdpsol_sedumi.PackagesAvailable
                message = sprintf('The following required packages are missing: %s.', strjoin(test_sdpsol_sedumi.MissingPackages, ', '));
                testCase.assumeTrue(test_sdpsol_sedumi.PackagesAvailable, message);
            end
        end
    end

    methods (TestParameterDefinition,Static)

       function [sdp, opts] = initializeTestData()

            % only run initialization if all required packages are available
            if ~test_sdpsol_sedumi.PackagesAvailable
                sdp  = {[]};
                opts = {[]};
                return;
            end

            % set seed
            rng(0)
            
            % basic data
            A = [-1 2 0;-3 -4 1;0 0 -2];
            x = casadi.SX.sym('x', 3, 3);
            
            % define f,g,x
            sdp = {struct('f',x(1,1),'g',[trace(x); vec(-A'*x-x*A)],'x',x(:))};

            % define cones
            opts = {struct('Kx',struct('psd',3),'Kc',struct('lin',1,'psd',3))};

       end

    end


    methods (Test)

        function solve_sdp(testCase,sdp,opts)
            
            % initialize solver
            S = casos.sdpsol('S','sedumi',sdp,opts);
            
            % solve with equality constraints
            sol = S('lbg', 1, 'ubg',1);

            if S.stats.UNIFIED_RETURN_STATUS == "SOLVER_RET_SUCCESS"
                actSolution = 1;
                refSolution = 1;
            else
                refSolution = 1;
                actSolution = inf;
            end

            % Perform assertions if needed
            testCase.verifyEqual(actSolution, refSolution ,"AbsTol",1e-12);

        end

    end
   
end