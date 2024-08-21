function [alpha_min] = computeMin_stepLength(obj,xi_k,p0,d_star,curr_conVio)

s_phi                = obj.opts.s_phi;
s_theta              = obj.opts.s_theta;
gamma_phi            = obj.opts.gamma_phi;
delta                = obj.opts.delta;

% compute alpha min
dirVecDeriv = full(casadi.DM(full(obj.nabla_xi_f(xi_k,p0))*d_star));
if dirVecDeriv < 0
    
    alpha_min = 0.05*min([ obj.opts.alpha_min, ( gamma_phi*curr_conVio )/abs(dirVecDeriv ), (delta*curr_conVio^s_theta )/(abs(dirVecDeriv)^s_phi) ]);
    
else
    
    alpha_min = obj.opts.alpha_min; % same value as gamma_theta in filter
    
end
end