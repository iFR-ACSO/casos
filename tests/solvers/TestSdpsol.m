% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

classdef TestSdpsol < TestSolver
% Test SDP solver interface.

properties (Constant)
    ALL_SOLVERS = {'sedumi' 'mosek' 'scs' 'clarabel'};
    EXCLUDED_SOLVERS = {'clarabel', 'scs'};
end

properties (TestParameter)
    solver = setdiff(TestSdpsol.ALL_SOLVERS, TestSdpsol.EXCLUDED_SOLVERS);
end

methods (Test)
    function test_sdp(test_case, solver)
        % Solve simple SDP with only linear and psd cones with sedumi
        % see https://yalmip.github.io/tutorial/semidefiniteprogramming/
        %               
        % Test problem: 
        %
        %   min       f         |   f   = x(1)         
        %   s.t.    x in K_x    |   K_x = single psd cone
        %          Ax in K_c    |   K_c = single psd cone 
        %                               
        test_case.check_required_package(solver);

        % basic data
        A = [-1 2 0;-3 -4 1;0 0 -2];
        x = casadi.SX.sym('x', 3, 3);
        
        % define f,g,x
        sdp = struct('f',x(1,1),'g',[trace(x); vec(-A'*x-x*A)],'x',x(:));

        % define cones
        opts = struct('Kx',struct('psd',3),'Kc',struct('lin',1,'psd',3));

        % initialize solver
        S = casos.sdpsol('S',solver,sdp,opts);
        
        % solve with equality constraints
        sol = S('lbg', 1, 'ubg',1);

        % perform assertions
        test_case.verifyEqual(S.stats.UNIFIED_RETURN_STATUS,casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS);
        test_case.verifyEqual(full(sol.g(1)),1,"AbsTol",1e-6);
        test_case.verifyGreaterThanOrEqual(eig(full(reshape(sol.x,3,3))),0);
        test_case.verifyGreaterThanOrEqual(eig(full(reshape(sol.g(2:10),3,3))),0);
    end

    function test_socp(test_case, solver)
        % Solve simple second order cone programm with mosek
        %               
        % Test problem: 
        %
        %   min       f
        %   s.t.    x in K_x    |   K_x = linear cone
        %          Ax in K_c    |   K_c = two socp cones
        %                               
        test_case.check_required_package(solver);

        % basic data
        t = (0:0.02:2*pi)';
        x = (1:6)';
        A = sin(t * x');
        y = A*x;
        
        x = casadi.SX.sym('x', 6, 1);
        u = casadi.SX.sym('u', 1, 1);
        v = casadi.SX.sym('v', 1, 1);
        
        % define f,g,x
        sdp = struct('f', u+v,'g',[u; y-A*x; v; x],'x', [x; u; v]);

        % define cones
        opts = struct('Kx',struct('lin',8),...
                      'Kc',struct('lor',[length(y)+1; length(x)+1]));

        % initialize solver
        S = casos.sdpsol('S',solver,sdp,opts);
        
        % solve without variable bounds
        sol = S('lbx',-inf,'ubx',inf);

        % separate solution
        G = mat2cell(sol.g,opts.Kc.lor);

        % perform assertions
        test_case.verifyEqual(S.stats.UNIFIED_RETURN_STATUS,casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS);
        test_case.verifyGreaterThanOrEqual(full(G{1}(1)),norm(full(G{1}(2:end)),2) - 1e-5);
        test_case.verifyGreaterThanOrEqual(full(G{2}(1)),norm(full(G{2}(2:end)),2) - 1e-5);
    end
end

end
