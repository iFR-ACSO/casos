%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Short Description: 
% 
% Call the actual quadratic subproblem (SDP with quadratic cost) of the
% form 
%   min 1/2*d'*B*d + ∇f'*d
%   st.  ∇g(x)⋅d + g(x+d0) - ∇g(x)⋅d0 in Σ[x] 
%       x in R[x] 
%
% where d = (x^* - x_k) is the search direction, d0 is the search direction
% (coming from the actual Q-SDP) i.e., searach direction we want to correct,
% B is the approximation of the Hessian, ∇(⋅) the first-derivative 
% and (⋅)^* is the optimal solution.
%
% Remark: All decision variables are setup as regular polynomials. SOS
% polynomials are enforced via constraints!
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_star_soc,dual_star_soc,skip_soc] = solve_Q_SDP_soc(obj,x_k,x_star,p0,Bk,args)

import casos.package.UnifiedReturnStatus


% pre-pare input for Q-SDP for soc
args{2}  = [p0; x_k; Bk(:);x_star];
sol_soc   = eval_on_basis(obj.solver_soc, args);


%% check solution status
switch (obj.solver_soc.get_stats.UNIFIED_RETURN_STATUS)

    case UnifiedReturnStatus.SOLVER_RET_SUCCESS    % optimal solution found
        x_star_soc    = sol_soc{1};
        dual_star_soc    = sol_soc{3};
        skip_soc      = 0;

    case UnifiedReturnStatus.SOLVER_RET_UNKNOWN    % feasible solution but not optimal

        x_star_soc    = sol_soc{1};
          dual_star_soc    = sol_soc{3};
        skip_soc      = 0;


    otherwise

        x_star_soc  = sol_soc{1};
          dual_star_soc    = sol_soc{3};
        skip_soc    = 1;
        
end

end