% ========================================================================
%
% Test Name: test_sdpsol.m
%
% Test Description: 
%   > Check for a simple SDP with only linear and psd cones if the
%   interface with mosek and sedumi works.
%               
% Test Procedure: ToDo
%                               
% Date: 09/06/2024 
%
% ========================================================================

classdef test_sdpsol < matlab.unittest.TestCase

    properties (TestParameter)
        sdp  
        opts 
    end

    properties (Access = private, Constant)
        PackagesAvailable   = test_sdpsol.checkRequiredPackages(1);
        MissingPackages     = test_sdpsol.checkRequiredPackages(2);
    end
    
    methods (Static)
        function output = checkRequiredPackages(out_select)
            % Define a list of required packages and their check functions
            packages = {'sedumi', 'multipoly', 'sosopt', 'mosek'};
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
            if ~test_sdpsol.PackagesAvailable
                message = sprintf('The following required packages are missing: %s.', strjoin(test_sdpsol.MissingPackages, ', '));
                testCase.assumeTrue(test_sdpsol.PackagesAvailable, message);
            end
        end
    end

    methods (TestParameterDefinition,Static)

       function [sdp, opts] = initializeTestData()

            % only run initialization if all required packages are available
            if ~test_sdpsol.PackagesAvailable
                sdp  = {[]};
                opts = {[]};
                return;
            end

            % set seed
            rng(0)
        
            % dimension of SDP
            ndim = 2;
        
            % decision variables
            x = casadi.SX.sym('x', 4, 1);
        
            % create random SDP data
            A = randn(ndim, ndim);
            B = randn(ndim, ndim);
            A = A+A';
            B = B+B';
        
            % define SDP problem
            sdp.f = x(1);
            sdp.g = eye(ndim) + x(1)*A + x(2)*B;
            sdp.g = [sdp.g(:)];
            sdp.x = x;
       
            % define cones
            opts.Kx = struct('lin', 4);
            opts.Kc = struct('psd', ndim);

            sdp  = {sdp};
            opts = {opts};
       end

    end


    methods (Test)

        function test_sdpsol_mosek(testCase,sdp,opts)
           
            % initialize solver
            S = casos.sdpsol('S','mosek',sdp,opts);
            
            % solve with equality constraints
            sol = S('lbx',-inf,'ubx',+inf);

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

        function test_sdpsol_sedumi(testCase,sdp,opts)
            
            % initialize solver
            S = casos.sdpsol('S','sedumi',sdp,opts);
            
            % solve with equality constraints
            sol = S('lbx',-inf,'ubx',+inf);

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