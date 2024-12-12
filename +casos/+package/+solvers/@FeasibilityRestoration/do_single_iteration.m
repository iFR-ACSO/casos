function [sol_iter,sol_qp,feas_res_flag,info,obj,filter,Bk] = do_single_iteration(obj, ...
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

%% evaluate convex SOS problem and check feasibility
[x_star,dual_star,sol_iter,sol_qp,feas_res_flag,info] = solve_Q_SDP(obj,iter,x_k,p0,Bk,args,info);


if feas_res_flag
    return % leave single iteration and invoke feasibility restoration
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

theta_min = min(filter(1,1),1)*1e-4;

while true
    
    % compute search direction
    dk = (x_star    - x_k);
    
    % check filter acceptance
    [x_k1,theta_x_k1,f_x_k1 ,filter_Acceptance] = check_filter_acceptance(obj,filter,alpha,x_k,dk,p0,args);
    
    if filter_Acceptance % acceptable to filter

         % check sufficient decrease
        [suffDecrease_flag,f_type,amijo] = chechSuffDecrease(obj,alpha,x_star,x_k,p0,theta_xk,theta_x_k1, f_x_k1,f_xk, ...
                                                             theta_min,eta,delta,gamma_theta,gamma_phi,s_phi,s_theta);
    
        
        if suffDecrease_flag
            % accept and leave while loop
            break
        end
        
        % invoke soc to avoid maratos effect, if necessary
        if ~suffDecrease_flag && alpha == 1 && theta_x_k1 >= theta_xk

           % compute correct search direction, compute new trial point and
           % check for filter acceptance
           x_k1 = second_order_correction(obj,x_k,x_star,p0,Bk,args,filter,alpha,theta_xk,f_xk,...
                                      theta_min,eta,delta,gamma_theta,gamma_phi,s_phi,s_theta);
           
           % leave while loop if corrected step is acceptable to filter
           if ~isempty(x_k1)
               break
           end
              
        end


    else % not acceptable 
        
        % invoke soc to avoid maratos effect, if necessary
        if ~filter_Acceptance && alpha == 1 && theta_x_k1 >= theta_xk
            
           % compute correct search direction, compute new trial point and
           % check for filter acceptance
           x_k1 = second_order_correction(obj,x_k,x_star,p0,Bk,args,filter,alpha,theta_xk,f_xk,...
                                          theta_min,eta,delta,gamma_theta,gamma_phi,s_phi,s_theta);

           if ~isempty(x_k1)
               break
           end
              
        end

    end
    
    % if not acceptable to filter and no sufficient progress update alpha
    alpha = 1/2*alpha;
    
    % compute alpha_min

    if alpha < alpha_min
       % invoke feasibility restoration if step-length is below minimum
       feas_res_flag = 2;
       break
    end
        

end % end of while-loop

% augment filter if not in forbidden region
if filter_Acceptance 

   % only augment if either f-type or amijo are not fulfilled
   if ~f_type || ~amijo

     % augment filter with a small margin
     filter = [filter;[theta_xk*(1-gamma_theta), f_xk - gamma_phi*theta_xk]];

   end

end

% update dual variables
dual_k1 = full(dual_k + alpha*(dual_star - dual_k));

% output of current iterate
sol_iter.x_k1       = x_k1;
sol_iter.dual_k1    = dual_k1;
sol_iter.theta_x_k1 = theta_x_k1;
sol_iter.f_x_k1     = f_x_k1;
sol_iter.alpha_k    = alpha;


Bk = damped_BFGS(obj,Bk,x_k,p0,sol_iter);

end


