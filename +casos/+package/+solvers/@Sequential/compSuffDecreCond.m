function [accept,Filter] = compSuffDecreCond(obj,p0,xi_k1,xi_k,alpha_k,old_conVio,predictedConVio,predictedcost ,s_phi,s_theta,eta,gamma_phi,gamma_theta,Filter)
                
                % values that are needed several times; compute only once
                searchDir  = obj.delta_search(xi_k1,xi_k);
                nablaF     = obj.nabla_f(xi_k,p0);
                theta_min  = min(1e-4,Filter(1,2));
                delta = 1;

                % sufficient reduction in cost and constraint violation
                suff_dec_conVio = old_conVio < theta_min;
                suff_dec_cost   = double(nablaF'*searchDir) < 0 && double( alpha_k*(- nablaF'*searchDir )^s_phi ) > delta*old_conVio^s_theta;
                    

                % check for f-type switching
                if suff_dec_conVio && suff_dec_cost
                    
                    % check if Amijo-condition is fulfilled
                    if  double(obj.cost(xi_k1,p0)) <= double(obj.cost(xi_k,p0)) + double(alpha_k*eta*nablaF'*searchDir)
                         accept = 1;
                   
                    else
                        
                        accept = 0;
                    end

                else  % if descent condition is not fulfilled, we know that the current iterate is acceptable to the filter
                    Filter       = vertcat(Filter, [predictedcost, predictedConVio]);
                    % Filter       = vertcat(Filter, [double(obj.cost(xi_k,p0)), old_conVio ]);

                    % sufficient reduction in constraint violation or in cost
                    if predictedConVio <= (1-gamma_theta)*old_conVio || predictedcost  <= double(obj.cost(xi_k,p0)-gamma_phi*old_conVio)
                        
                          accept = 1;
                    else
                          accept = 0;
                    end

                end


end