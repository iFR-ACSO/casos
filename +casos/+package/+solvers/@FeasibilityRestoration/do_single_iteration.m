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
alpha_max   = obj.opts.filter_struct.alpha_max;

gamma_theta = obj.opts.filter_struct.gamma_theta ;
gamma_phi   = obj.opts.filter_struct.gamma_phi ;

LangrangeFilter = obj.opts.filter_struct.LangrangeFilter;


alpha = alpha_max;

% Backtracking linesearch
while true
    
    % compute search direction
    dk   = (x_star    - x_k);
    dkl  = (dual_star-dual_k);
     
    % check, if search direction becomes too small
    if max(abs(full(dk))./(1+abs(full(x_k)))) < 10*eps
        stop = 1;
    end

    % check filter acceptance
    [x_k1,theta_x_k1,f_x_k1 ,filter_Acceptance,sol_convio] = check_filter_acceptance(obj,filter,alpha,x_k,dual_k,dk,dkl,p0,args);
    
     dual_k1 = dual_k + alpha*(dual_star-dual_k);
    
     if LangrangeFilter
         L_k1= full(obj.L(x_k1,p0,dual_k1));
         L_k = full(obj.L(x_k,p0,dual_k));
     else
         L_k1 = [];
         L_k  = [];
     end

    if filter_Acceptance 

        % check sufficient decrease
        [suffDecrease_flag,f_type,amijo,filter] = checkSuffDecrease(obj,alpha,x_star,x_k,p0,theta_xk,theta_x_k1, f_x_k1,f_xk,L_k1,L_k,dual_k,filter);
    
        if suffDecrease_flag
            % accept and leave while loop
            break
        end
        
        % invoke soc to avoid maratos effect, if necessary
        if ~suffDecrease_flag && alpha == 1  && theta_x_k1 >= theta_xk

           % compute corrected search direction, compute new trial point and check for filter acceptance
           [x_k1,suffDecrease_flag,f_type,amijo] = second_order_correction(obj,x_k,x_star,p0,Bk,args,filter,alpha,theta_xk,f_xk,dkl,L_k1,L_k,dual_k);

           % leave while loop if corrected step is acceptable to filter
           if suffDecrease_flag
               break % means we found a solution with SOC; leave loop
           end

        end


    else % not acceptable to filter
        
        % invoke soc to avoid maratos effect, if necessary
        if ~filter_Acceptance && alpha == 1 && theta_x_k1 >= theta_xk
            
           % compute corrected search direction, compute new trial point and check for filter acceptance
           [x_k1,suffDecrease_flag,f_type,amijo] = second_order_correction(obj,x_k,x_star,p0,Bk,args,filter,alpha,theta_xk,f_xk,dkl,L_k1,L_k,dual_k);

            if suffDecrease_flag
               break % means we found a solution with SOC; leave loop
           end
              
        end

    end
    
    % if not acceptable to filter and/or no sufficient progress update alpha
    alpha = 1/2*alpha;
    
    % compute alpha_min
    alpha_min = compute_alpha_min(obj,x_k,dk,p0,theta_xk,filter);
    
    if alpha < alpha_min  || alpha <= eps
       % invoke feasibility restoration if step-length is below minimum
       feas_res_flag = 2;
       break
    end
        

end % end of backtracking linesearch


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


% output of current iterate
sol_iter.x_k1       = x_k1;
sol_iter.dual_k1    = dual_k1;
sol_iter.theta_x_k1 = theta_x_k1;
sol_iter.f_x_k1     = f_x_k1;
sol_iter.alpha_k    = alpha;
sol_iter.dual_qp    = dual_star;
info{iter}.constraint_violation = sol_convio;
% update quasi-Newton/exact Hessian
% if strcmp(obj.opts.hessian_approx,'BFGS')
        % damped BFGS
    Bk = damped_BFGS(obj,Bk,x_k,p0,sol_iter,iter);
% elseif strcmp(obj.opts.hessian_approx,'Levenberg')
%     % see Betts book
%     H  = full(obj.hess_fun(x_k1,p0,dual_k1));
%     Bk = casos.package.solvers.SequentialCommon.hessian_regularization(H);
% elseif strcmp(obj.opts.hessian_approx,'EigValReg')
%     % own regularization method
%     H = full(obj.hess_fun(x_k1,p0,dual_k1));
%     Bk = casos.package.solvers.SequentialCommon.regularize_Hessian(H);
% elseif strcmp(obj.opts.hessian_approx,'Mirroring')
%     % Verschueren
%     H = full(obj.hess_fun(x_k1,p0,dual_k1));
%     [L,D] = ldl(H);
% 
%     D(D<0) = 1e-6;
%     Bk = L*abs(D)*L';
% end


end
