function [accept] = compSuffDecreCond_SOC(obj,p0,x_soc_trial,xi_k1,xi_k,alpha_k,alpha_k0,old_conVio,predictedConVio,predictedcost ,s_phi,s_theta,eta,gamma_phi,gamma_theta,Filter)
                
                % IMPORTANT NOTE: An additional suff. decrease condition file was generated
                % due to the explanation given in [Wächter] below eq. 28
                % the current
            
                % values that are needed several times; compute only once
                searchDir  = obj.delta_search(x_soc_trial,xi_k);
                nablaF     = obj.nabla_f(xi_k,p0);
                theta_min  = 1e-4*max(1,Filter(1,2));
                delta = 1;

                % sufficient reduction in cost and violation
                suff_dec_conVio = old_conVio <= theta_min;

                % use alpha_k0 in suff. decrease in cost; use original
                % search direction as explained by Wächter
                suff_dec_cost   = (double(nablaF'*searchDir) < 0 && double( alpha_k*(- nablaF'*searchDir )^s_phi ) > delta*old_conVio^s_theta);

                if suff_dec_conVio && suff_dec_cost
                    
                    % check if Amijo-condition is fulfilled; ; use original
                    % search direction as explained by Wächter
                    % only on the left-hand side we use trial from SOC; not
                    % sure about which alpha should be used on the right
                    % hand side; currently alpha_soc
                    if  double(obj.cost(x_soc_trial,p0)) <= double(obj.cost(xi_k,p0)) + double(alpha_k*eta*nablaF'*searchDir)
                         accept = 1;
                    else
                        accept = 0;
                    end

                elseif ~suff_dec_conVio || ~suff_dec_cost % if stronger condition can not be fulfilled check if sufficient decrease w.r.t. Filter

                    % sufficient reduction in constraint violation or in
                    % cost
                    if predictedConVio <= (1-gamma_theta)*old_conVio || predictedcost  <= double(obj.cost(xi_k,p0)-gamma_phi*old_conVio)
                          accept = 1;
                    else
                          accept = 0;
                    end
                    
                end


end