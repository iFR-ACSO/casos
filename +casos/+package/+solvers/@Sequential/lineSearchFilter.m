function [xi_k1, dual_k1, Filter, Linesearch_stall,alpha_k,predictedcost,predictedConVio]   = lineSearchFilter(obj,...
                                                                                                                     p0,...
                                                                                                                     Filter,...
                                                                                                                     xi_k,...
                                                                                                                     xi_plus,...
                                                                                                                     dual_k,...
                                                                                                                     dual_plus,...
                                                                                                                     old_conVio,...
                                                                                                                     H,...
                                                                                                                     args_solver)   

        %% backtracking line search filter (ipopt); see Biegler Algorithm 5.4
        alpha_max            = 1;
        alpha_min            = 0.0001;
        tau                  = 0.5;
        
        gamma_theta = 0.01;
        gamma_phi  =  1;

        Linesearch_stall = 0;
        
        % parameter from IPOPT default parameters: https://coin-or.github.io/Ipopt/OPTIONS.html
        s_phi   = 2.3;
        s_theta = 1.1;
        eta     = 1e-4;

        alpha_k = alpha_max;
		
        while alpha_k >= alpha_min

            % compute new trial point
            xi_k1   = xi_k   + alpha_k*(xi_plus-xi_k);
            
  
            % check if current iterate is acceptable to filter
            [filter_accept,Filter,predictedConVio,predictedcost] = checkFilterAcceptance(obj,Filter,xi_k1,p0,gamma_theta,gamma_phi);
           
            % if acceptable, then check if sufficient decrease conditions are fulfilled
            if filter_accept
                     
                    [accept,Filter] = compSuffDecreCond(obj,...
                                                         p0,...
                                                         xi_k1,...
                                                         xi_k,...
                                                         alpha_k,...
                                                         old_conVio,...
                                                         predictedConVio,...
                                                         predictedcost,...
                                                         s_phi,...
                                                         s_theta,...
                                                         eta,...
                                                         gamma_phi,...
                                                         gamma_theta,...
                                                         Filter);

          
                % inner-loop i.e. second-order-correction (SOC)
                if ~accept  
                    
                    [xi_k1,alpha_k] = computeSOC(obj,args_solver,p0,xi_k,xi_plus,H,s_phi,s_theta,eta,gamma_phi,gamma_theta,Filter,old_conVio,tau,alpha_k);

                    if ~isempty(xi_k1)
                        break
                    end

                else

                      
                    break
                end
                
            else
                % if not acceptable to filter, check SOC
                if alpha_k == 1 
                    [xi_k1,alpha_k] = computeSOC(obj,args_solver,p0,xi_k,xi_plus,H,s_phi,s_theta,eta,gamma_phi,gamma_theta,Filter,old_conVio,tau,alpha_k);
                      
                    if ~isempty(xi_k1)
                        break
                    end

                else
                     % not acceptable to filter with lineserach
                    alpha_k = alpha_k*tau ;
                end
           
            end
                
               
        end % end-while outer backtracking
        
 
         if alpha_k < alpha_min
             %% go to feasibility restoration phase
             % x0 must contain current infesible initial guess plus an
             % initial guess for the "projection" multiplier
             % % lower and upper bounds needed
             % sol_feas =    obj.FeasResPhase('x0',x0, ...
             %                               'lbx',[Vlb;s2lb;slb], ...
             %                                'ubx',[Vub;s2ub;sub]);

             Linesearch_stall = 0;


         end
        
         % compute Langrange variables
         if ~isempty(dual_k)
                dual_k1 = dual_k + alpha_k*(dual_plus-dual_k);
         else
                dual_k1 = dual_plus;
         end

end
