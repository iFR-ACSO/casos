%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Short Description: 
% 
% Call the actual quadratic subproblem (SDP with quadratic cost) of the
% form 
%   min 1/2*d'*B*d + ∇f'*d
%   st. g(x) + ∇g(x)⋅d in Σ[x] 
%       x in R[x] 
%
% where d = (x^* - x_k) is the search direction, B is the approximation of 
% the Hessian, ∇(⋅) the first-derivative and (⋅)^* is the optimal solution.
%
% Remark: All decision variables are setup as regular polynomials. SOS
% polynomials are enforced via constraints!
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_star,dual_star,sol_iter,sol_qp,feas_res_flag,info] = solve_Q_SDP(obj, ...
                                                                             iter, ...
                                                                             x_k, ...
                                                                             p0, ...
                                                                             Bk, ...
                                                                             args, ...
                                                                             info)

import casos.package.UnifiedReturnStatus

% pre-pare input for Q-SDP
args{2}  = [p0; x_k; Bk(:)];
% adjust bounds
% args{3}  = args{3} - x_k;
% args{4}  = args{4} - x_k;
measSolver =tic;
sol_qp   = eval_on_basis(obj.solver_convex, args);
solveTimeSDP = toc(measSolver);
% store iteration info
info{iter} = obj.solver_convex.get_stats  ;
info{iter}.timeStats.mainSDP = solveTimeSDP;

%% check solution status
switch (obj.solver_convex.get_stats.UNIFIED_RETURN_STATUS)
    
     % optimal solution found (or at least almost feasible)
    case UnifiedReturnStatus.SOLVER_RET_SUCCESS   
        
        x_star    = sol_qp{1};
        dual_star = sol_qp{5};
        feas_res_flag = 0;
        sol_iter = [];
    
      % case where the solver could not decide or numerical issues
     case UnifiedReturnStatus.SOLVER_RET_LIMITED

         feas_res_flag = 1;
         x_star        = x_k;
         dual_star     = [];
         sol_iter      = [];

    % strictly infeasible 
    otherwise
        feas_res_flag = 1;
        x_star    = x_k;
        dual_star = [];
        sol_iter = [];

end

end