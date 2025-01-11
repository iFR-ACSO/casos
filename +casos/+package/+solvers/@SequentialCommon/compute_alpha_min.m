function alpha_min = compute_alpha_min(obj,x_k,dk,p0,s_phi,delta,theta_xk,s_theta,gamma_theta,gamma_phi,gamma_alpha,theta_min)

% descent direction
nabla_f_dir = full(obj.eval_gradCost(x_k,p0)*dk);

% theta_xk can be zero. in that case
if nabla_f_dir  < 0 && theta_xk <= theta_min

   alpha_min = min([gamma_theta, (gamma_phi*theta_xk)/(-nabla_f_dir), (delta*theta_xk^s_theta)/((-nabla_f_dir)^s_phi) ]);

elseif nabla_f_dir  < 0 && theta_xk > theta_min
    
    alpha_min =  min([gamma_theta,(gamma_phi*theta_xk)/(-nabla_f_dir)]);

else
    
   alpha_min = gamma_theta;

end

 alpha_min = gamma_alpha*alpha_min;

end
