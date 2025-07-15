%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Short Description: 
% 
% Compute minimum step-length based on sufficient descrease condition, i.e.
% check if these can be fulfilled (based on current iterate) if alpha is
% small enough. For a derivation see equation (23) in
%
% Wächter, A. and Biegler, L. - On the implementation of an interior-point 
% filter line-search algorithm for large-scale nonlinear programming,  
% Mathematical Programming, 2006,doi: 10.1007/s10107-004-0559-y
% 
% Note: alpha_min can be zero!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alpha_min = compute_alpha_min(obj,x_k,dk,p0,theta_xk,filter)

% extract options from struct; for readability oc code
s_theta     = obj.opts.filter_struct.s_theta ;
s_phi       = obj.opts.filter_struct.s_phi ;
gamma_theta = obj.opts.filter_struct.gamma_theta ;
gamma_phi   = obj.opts.filter_struct.gamma_phi ;

gamma_alpha     = obj.opts.filter_struct.gamma_alpha;

delta           = obj.opts.filter_struct.delta;

LangrangeFilter = obj.opts.filter_struct.LangrangeFilter;

theta_min = min(filter(1,1),1)*obj.opts.filter_struct.theta_min;

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

% scale alpha min.
alpha_min = gamma_alpha*alpha_min;

end
