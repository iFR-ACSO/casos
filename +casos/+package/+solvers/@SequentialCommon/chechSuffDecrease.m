function [suffDecrease_flag,f_type,amijo,filter] = chechSuffDecrease(obj,alpha,x_star,x_k,p0,theta_xk,theta_x_k1, f_x_k1,f_xk, ...
                                                              theta_min,eta,delta,gamma_theta,gamma_phi,s_phi,s_theta,filter)


% search direction
dk = x_star    - x_k;

% descent direction
nabla_f_dir = full(obj.eval_gradCost(x_k,p0)*dk);

% f-type switching: descent direction + progress in cost better than con.
% violation
f_type = nabla_f_dir  < 0 && alpha*(-nabla_f_dir)^s_phi > delta*theta_xk^s_theta;

% initialize amijo boolean
amijo  = 0;

% only augment filter if not f-type
if ~f_type
     % augment filter with a small margin
     filter = [filter;[theta_xk*(1-gamma_theta), f_xk - gamma_phi*theta_xk]];
end

if theta_xk < theta_min && f_type
    % check amijo-condition
    amijo = f_x_k1 <= f_xk + eta*alpha*nabla_f_dir;

    if amijo
        suffDecrease_flag = 1;
    else
        suffDecrease_flag = 0;
    end

else
  
    % check progress w.r.t. previous iteration
    if theta_x_k1 <= (1-gamma_theta)*theta_xk || f_x_k1 <= f_xk - gamma_phi*theta_xk
        suffDecrease_flag = 1;
    else
        suffDecrease_flag = 0;
    end % check if soc shall be used or if alpha needs adjustment

end % end of check progress w.r.t. current iterate

end % end of function