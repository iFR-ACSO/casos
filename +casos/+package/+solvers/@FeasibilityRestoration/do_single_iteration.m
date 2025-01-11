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
s_theta     = 0.9;
s_phi       = 2;
gamma_theta = 1e-5;
gamma_phi   = 1e-5;

gamma_alpha = 0.05;

delta       = 1;
eta         = 1e-4;
LangrangeFilter = 0; 
theta_min = min(filter(1,1),1)*1e-4;

% Backtracking linesearch
while true
    
    % compute search direction
    dk   = (x_star    - x_k);
    dkl  = (dual_star-dual_k);
     
    % check filter acceptance
    [x_k1,theta_x_k1,f_x_k1 ,filter_Acceptance] = check_filter_acceptance(obj,filter,alpha,x_k,dual_k,dk,dkl,p0,args);
    
     dual_k1 = dual_k + alpha*(dual_star-dual_k);
    
     if LangrangeFilter
         L_k1= full(obj.L(x_k1,p0,dual_k1));
         L_k = full(obj.L(x_k,p0,dual_k));
     else
         L_k1 = [];
         L_k  = [];
     end

    if filter_Acceptance  %&& Wolfe % acceptable to filter

        % check sufficient decrease
        [suffDecrease_flag,f_type,amijo,filter] = chechSuffDecrease(obj,alpha,x_star,x_k,p0,theta_xk,theta_x_k1, f_x_k1,f_xk,L_k1,L_k,dual_k, ...
                                                                    theta_min,eta,delta,gamma_theta,gamma_phi,s_phi,s_theta,filter);
    
        
        % either sufficient decrease in cost or constraint violation
        if suffDecrease_flag
            % accept and leave while loop
            break
        end
        
        % invoke soc to avoid maratos effect, if necessary
        if ~suffDecrease_flag && alpha == 1 

           % compute corrected search direction, compute new trial point and check for filter acceptance
           [x_k1,suffDecrease_flag,f_type,amijo] = second_order_correction(obj,x_k,x_star,p0,Bk,args,filter,alpha,theta_xk,f_xk,dkl,L_k1,L_k,dual_k,...
                                                                            theta_min,eta,delta,gamma_theta,gamma_phi,s_phi,s_theta);

           % leave while loop if corrected step is acceptable to filter
           if ~isempty(x_k1)
               break % means we found a solution with SOC; leave loop
           end

        end


    else % not acceptable to filter
        
        % invoke soc to avoid maratos effect, if necessary
        if ~filter_Acceptance && alpha == 1 %&& theta_x_k1 >= theta_xk
            
           % compute corrected search direction, compute new trial point and check for filter acceptance
           [x_k1,suffDecrease_flag,f_type,amijo] = second_order_correction(obj,x_k,x_star,p0,Bk,args,filter,alpha,theta_xk,f_xk,dkl,L_k1,L_k,dual_k,...
                                                                             theta_min,eta,delta,gamma_theta,gamma_phi,s_phi,s_theta);

           if ~isempty(x_k1)
               break % means we found a solution with SOC; leave loop
           end
              
        end

    end
    
    % if not acceptable to filter and/or no sufficient progress update alpha
    alpha = 1/2*alpha;
    
    % compute alpha_min
    alpha_min = compute_alpha_min(obj,x_k,dk,p0,s_phi,delta,theta_xk,s_theta,gamma_theta,gamma_phi,gamma_alpha);
    
    if alpha < alpha_min 
       % invoke feasibility restoration if step-length is below minimum
       feas_res_flag = 2;
       break
    end
        

end % end of while-loop


% only augment if new iterate is accepted and either f-type or amijo are not fulfilled
if suffDecrease_flag && (~f_type || ~amijo)

    if LangrangeFilter
        % augment filter with a small margin
        filter = [filter;[theta_xk*(1-gamma_theta), L_k - gamma_phi*theta_xk]];
    else
        % augment filter with a small margin
        filter = [filter;[theta_xk*(1-gamma_theta), f_xk - gamma_phi*theta_xk]];
    end


end


% update dual variables
% dual_k1 = full(dual_k + alpha*(dual_star - dual_k));

% dual_k1 = dual_star;
% output of current iterate
sol_iter.x_k1       = x_k1;
sol_iter.dual_k1    = dual_k1;
sol_iter.theta_x_k1 = theta_x_k1;
sol_iter.f_x_k1     = f_x_k1;
sol_iter.alpha_k    = alpha;
sol_iter.dual_qp    = dual_star;





% Bk = (H0+H0')/2;

% update BFGS matrix for next iteration
% % if colFlag
    Bk = damped_BFGS(obj,Bk,x_k,p0,sol_iter,iter);
% else
    % Bk = H;
% end

end
