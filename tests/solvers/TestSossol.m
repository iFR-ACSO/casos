% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

classdef TestSossol < TestSolver
% Test sum-of-squares solver interface.

properties (TestParameter)
    solver = {'sedumi' 'mosek' 'scs' 'clarabel'};
end

methods (Test)
    function solve_sos(test_case, solver)
        %   > Solve simple sos program
        %               
        % Test Procedure: 
        %
        %   min       f
        %   s.t.    x in K_x    |   K_x = linear cone
        %          Ax in K_c    |   K_c = single sos cone
        %
        test_case.check_required_package(solver);

        % indeterminate variable
        x = casos.Indeterminates('x');

        % some polynomial
        f = x^4 + 10*x;

        % scalar decision variable
        g = casos.PS.sym('g');

        % define f,g,x
        sos = struct('x',g,'f',g,'g',f+g);

        % define cones
        opts = struct('Kc',struct('sos',1));

        % initialize solver
        S = casos.sossol('S',solver,sos,opts);
        
        % solve with equality constraints
        sol = S();

        % perform assertions
        test_case.verifyEqual(S.stats.UNIFIED_RETURN_STATUS,casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS);
        test_case.verifyEqual(full(sol.f),full(sol.x));
        test_case.verifyGreaterThanOrEqual(full(sol.f),10.179);
    end
end

end
