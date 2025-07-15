%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Short Description: 
% 
% Implementation of the damped Broyden–Fletcher–Goldfarb–Shanno (BFGS) method
% to estimate the Hessian. 
% For details see:
%    
% Nocedal, J., Wright, S. - Numerical Optimization, Springer, 2006, on page
% 536-537
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Bk = damped_BFGS(obj,Bk,x_k,p0,sol_iter,iter)

% update BFGS
s_val = full( obj.eval_s(x_k,sol_iter.x_k1,p0) );
y_val = full( obj.eval_y(x_k,sol_iter.x_k1,p0,sol_iter.dual_k1) );
        

 % Powell damping
 if s_val'*y_val >= 0.2*s_val'*Bk*s_val 
     theta = 1;
 else
     theta = (0.8*s_val'*Bk*s_val)/(s_val'*Bk*s_val-s_val'*y_val);
 end

 r_val = full( obj.eval_r(Bk,theta,y_val,s_val) );

 % update Bk
 % B = B + (r*r')/(s'*r) - (B*s*s'*B)/(s'*B*s);
 Bk =  full( obj.damped_BFGS(Bk,r_val,s_val) );
 
 % check symmetry of Hessian and ajdust if necessary
 Bknonsym = norm(Bk-Bk',inf);
 
 if Bknonsym > eps
     Bk = (Bk+Bk')/2;
 end
 
 % check if Hessian is PD, if not apply perturbation
 [~,cholFlag] = chol(Bk);

 if cholFlag
     dE=abs(eig(Bk,'nobalance'));
     pert=(min(dE)+2*eps*max(dE))/(1-2*eps);
    Bk=Bk + pert*eye(length(Bk));
 end
 
           
end
