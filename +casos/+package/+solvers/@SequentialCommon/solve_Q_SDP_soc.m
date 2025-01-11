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