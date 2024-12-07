function Bk = damped_BFGS(obj,Bk,x_k,p0,sol_iter)

       % update BFGS
       s_val = full( obj.eval_s(x_k,sol_iter.x_k1,p0) );
       y_val = full( obj.eval_y(x_k,sol_iter.x_k1,p0,sol_iter.dual_k1) );
        
        if s_val'*Bk*s_val < eps || s_val'*y_val < eps
            % skip
        else

            % Powell damping
            if s_val'*y_val >= 0.2*s_val'*Bk*s_val 
                theta = 1;
            else
                theta = (0.8*s_val'*Bk*s_val)/(s_val'*Bk*s_val-s_val'*y_val);
            end
           
            r_val = full( obj.eval_r(Bk,theta,y_val,s_val) );
            
            % B = B + (r*r')/(s'*r) - (B*s*s'*B)/(s'*B*s);
            Bk =  full( obj.damped_BFGS(Bk,r_val,s_val) );
            
            % check symmetry of Hessian and ajdust if necessary
            Bknonsym = norm(Bk-Bk',inf);

            if Bknonsym > eps
                Bk = (Bk+Bk')/2;
            end
        end
end
