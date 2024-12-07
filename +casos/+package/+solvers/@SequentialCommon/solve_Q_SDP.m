function [x_star,dual_star,sol_iter,sol_qp,feas_res_flag,info] = solve_Q_SDP(obj,iter,x_k,p0,Bk,args,info)

import casos.package.UnifiedReturnStatus

% pre-pare input for Q-SDP
args{2}  = [p0; x_k; Bk(:)];
sol_qp   = eval_on_basis(obj.solver_convex, args);

% store iteration info
info{iter} = obj.solver_convex.get_stats;

%% check solution status
switch (obj.solver_convex.get_stats.UNIFIED_RETURN_STATUS)

    case UnifiedReturnStatus.SOLVER_RET_SUCCESS    % optimal solution found
        x_star    = sol_qp{1};
        dual_star = sol_qp{5};
        feas_res_flag = 0;
        sol_iter = [];

    case UnifiedReturnStatus.SOLVER_RET_UNKNOWN    % feasible solution but not optimal
        x_star    = sol_qp{1};
        dual_star = sol_qp{5};
        feas_res_flag = 0;
        sol_iter = [];


    case UnifiedReturnStatus.SOLVER_RET_INFEASIBLE % actually this is otherwise-case
            feas_res_flag = 1;
        
        x_star = sol_qp{1};
        dual_star = sol_qp{5};
         % store already for next iteration
        sol_iter.x_k1      = [];
        sol_iter.dual_k1   = [];

        % compute constraint violation
        % work around to get correct arguments for constrained violation check
        args_conVio     =  args;
        args_conVio{2}  =  [p0; sol_qp{1}];
        args_conVio{3}  = -inf(obj.init_para.conVio.no_con,1);
        args_conVio{4}  =  inf(obj.init_para.conVio.no_con,1);
        
        % constraint violation at trial point
        sol_convio = eval_on_basis(obj.solver_conVio, args_conVio);
        
        theta_x_k1 = full(max(0,max(sol_convio{1})));
    
        % cost at trial point
        f_x_k1   = full(obj.eval_cost(sol_qp{1},p0));

        sol_iter.theta_x_k1 = theta_x_k1;
        sol_iter.f_x_k1     = f_x_k1;
    otherwise
        x_star    = sol_qp{1};
        dual_star = sol_qp{5};
        feas_res_flag = 0;
        sol_iter = [];
        


end

end