function [B] = dampedBFGS(obj,B,s,y)

% calculate damping parameter; Nocedal eq. (18.15)

if s'*B*s < 1e-14 || s'*y < 1e-14
    % skip
else
if full(casadi.DM( s'*y )) >= full(casadi.DM (0.2*s'*B*s ))
    theta = 1;
else
    theta = (0.8*s'*B*s)/(s'*B*s-s'*y);
end

% Powell damping
r = theta*y+(1-theta)*B*s;

% B = B + (r*r')/(s'*r) - (B*s*s'*B)/(s'*B*s);
B =  full(obj.BFGS_fun((B),(r),(s)));

if min(eig(B)) < 0
   warning('BFGS has negative eigenvalues!')

end

end
end