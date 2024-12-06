function [sol_iter,sol_qp,feas_res_flag,info,obj,filter] = do_single_iteration(obj, ...
                                                                              iter,...
                                                                              x_k,...
                                                                              dual_k,...
                                                                              theta_xk,...
                                                                              f_xk,...
                                                                              Bk,...
                                                                              p0,...
                                                                              args, ...
                                                                              filter, ...
                                                                              info)

import casos.package.UnifiedReturnStatus


%% evaluate convex SOS problem

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

    case UnifiedReturnStatus.SOLVER_RET_UNKNOWN    % feasible solution but not optimal
        x_star    = sol_qp{1};
        dual_star = sol_qp{5};
        feas_res_flag = 0;
    otherwise

    % case UnifiedReturnStatus.SOLVER_RET_INFEASIBLE % actually this is otherwise-case
        feas_res_flag = 1;

         % store already for next iteration
        sol_iter.x_k1      = sol_qp{1};
        sol_iter.dual_k1   = sol_qp{5};

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

        return
 
end


%% backtracking filter-linesearch

% To-Do parameter as options
alpha       = 1;
alpha_min   = 1e-6;
s_theta     = 1.1;
s_phi       = 2.3;
gamma_theta = 1e-3;
gamma_phi   = 1e-3;
delta       = 1;
eta         = 1e-4;

theta_min = max(filter(1,1),1)*1e-4;


while true

    % compute new solution candidate
    x_k1    = full(x_k     + alpha*(x_star    - x_k));
    
    % compute constraint violation
    % work around to get correct arguments for constrained violation check
    args_conVio     =  args;
    args_conVio{2}  =  [p0; x_k1];
    args_conVio{3}  = -inf(obj.init_para.conVio.no_con,1);
    args_conVio{4}  =  inf(obj.init_para.conVio.no_con,1);
    
    % constraint violation at trial point
    sol_convio = eval_on_basis(obj.solver_conVio, args_conVio);
    
    theta_x_k1 = full(max(0,max(sol_convio{1})));

    % cost at trial point
    f_x_k1   = full(obj.eval_cost(x_k1,p0));

    % check filter acceptance
    theta_l = filter(:,1);
    f_l     = filter(:,2);
    
    % new point lies in forbidden region if both are larger than filter entries
    dominance_bool      = [];
    dominance_bool(:,1) = theta_x_k1 >= theta_l;
    dominance_bool(:,2) = f_x_k1 >= f_l;
    
    
    if any(all(dominance_bool, 2)) % means not acceptable to filter

        if alpha == 1
            soc_flag = 1;
        else
         alpha = 1/2*alpha;

            if alpha < alpha_min
                feas_res_flag = 1;
                break
            end
        end

    else
        % check sufficient decrease

        % search direction
        dk = x_star    - x_k;
        
        % descent direction
        nabla_f_dir = full(obj.eval_gradCost(x_k,p0)*dk);

        % f-type switching
        f_type = nabla_f_dir  < 0 && alpha*(-nabla_f_dir)^s_phi > delta*theta_xk^s_theta;
        % initialize amigjo boolean
        amijo  = 0;

        if theta_xk < theta_min && f_type

            % check amijo-condition
            amijo = f_x_k1 <= f_xk + eta*alpha*nabla_f_dir;

            if amijo
                break % accept trial point (leave while-loop)
            else
                if alpha == 1
                    soc_flag = 1;
                else
                 alpha = 1/2*alpha;
        
                    if alpha < alpha_min
                        feas_res_flag = 1;
                        break
                    end
                end
            end
        
        else
            % check progress w.r.t. previous iteration
            if theta_x_k1 <= (1-gamma_theta)*theta_xk || f_x_k1 <= f_xk - gamma_phi*theta_xk
                break                 % accept trial point (leave while-loop)
            else
                % actually go to SOC
                if alpha == 1
                    soc_flag = 1;
                else

                    alpha = 1/2*alpha;
        
                    if alpha < alpha_min
                        feas_res_flag = 1;
                        break
                    end
                end % check if soc shall be used or if alpha needs adjustment
            end % end of check progress w.r.t. current iterate
        end 
    end % end of filter acceptance

    

end % end of while-loop

% augment filter if not in forbidden reagion
if ~any(all(dominance_bool, 2)) 

   % only augment if either f-type or amijo are not fulfilled
   if ~f_type || ~amijo

      % augment filter with a small margin
     filter = [filter;[theta_xk*(1-gamma_theta),f_xk - gamma_phi*theta_xk]];

   end
end


dual_k1 = full(dual_k + alpha*(dual_star - dual_k));

% output of current iterate
sol_iter.x_k1       = x_k1;
sol_iter.dual_k1    = dual_k1;
sol_iter.theta_x_k1 = theta_x_k1;
sol_iter.f_x_k1     = f_x_k1;
sol_iter.alpha_k    = alpha;

end


