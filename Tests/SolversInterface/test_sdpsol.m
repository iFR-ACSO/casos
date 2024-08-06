classdef test_sdpsol < matlab.unittest.TestCase

    properties (TestParameter)
        sdp  
        opts 
    end
    

    methods (TestParameterDefinition,Static)

       function [sdp, opts] = initializeTestData()

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