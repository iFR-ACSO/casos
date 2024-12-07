function  [x_k1,theta_x_k1,f_x_k1 ,filter_Acceptance] = check_filter_acceptance(obj,filter,alpha,x_star,x_k,p0,args)

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
        filter_Acceptance = 0;
    else
        filter_Acceptance = 1;
    end

end