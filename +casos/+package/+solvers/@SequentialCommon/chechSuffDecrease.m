function [suffDecrease_flag,f_type,amijo] = chechSuffDecrease(obj,alpha,x_star,x_k,p0,theta_xk,theta_x_k1, f_x_k1,f_xk, ...
                                                              theta_min,eta,delta,gamma_theta,gamma_phi,s_phi,s_theta)

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

end 