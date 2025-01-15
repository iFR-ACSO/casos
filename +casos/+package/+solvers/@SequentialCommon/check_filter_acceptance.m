function  [x_k1,theta_x_k1,f_x_k1 ,filter_Acceptance] = check_filter_acceptance(obj,filter,alpha,x_k,dual_k,dk,dkl,p0,args)
    
    % use Langrange function instead of cost function in filter
    LangrangeFilter = obj.opts.filter_struct.LangrangeFilter;

    % compute new solution candidate
    x_k1    = full(x_k     + alpha*dk);
    
    % compute constraint violation of new solution candidate
    args_conVio     =  args;
    args_conVio{2}  =  [p0; x_k1];
    args_conVio{3}  = -inf(obj.init_para.conVio.no_con,1);
    args_conVio{4}  =  inf(obj.init_para.conVio.no_con,1);
    
    % constraint violation at trial point
    sol_convio = eval_on_basis(obj.solver_conVio, args_conVio);
    
    theta_x_k1 = full(max(0,max(sol_convio{1})));

    % cost at trial point
    if LangrangeFilter 
          dual_k1 = full(dual_k  + alpha*dkl);
          L_k1     = full(obj.L(x_k1,p0,dual_k1));
    else
          f_x_k1   = full(obj.eval_cost(x_k1,p0));
    end
    
  

    %% check filter acceptance
    % get filter inputs
    theta_l = filter(:,1);
    f_l     = filter(:,2);
    
    % new point lies in forbidden region if both are larger than filter entries
    dominance_bool      = [];
    
    % boundary belongs to forbidden region i.e. >=
    dominance_bool(:,1) = theta_x_k1 >= theta_l; 
    
    if LangrangeFilter 
          dominance_bool(:,2) = L_k1     >= f_l;
    else
           dominance_bool(:,2) = f_x_k1  >= f_l;
    end
    % 
    
    % check pairs; if both one, means pair lies in forbidden region
    dominance_bool = all(dominance_bool, 2);
    
    if any(dominance_bool) % means not acceptable to filter
        filter_Acceptance = 0;
    else
        filter_Acceptance = 1;
    end

end