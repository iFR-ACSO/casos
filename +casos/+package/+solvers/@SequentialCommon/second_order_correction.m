function x_k1 = second_order_correction(obj,x_k,x_star,p0,Bk,args,filter,alpha,theta_xk,f_xk,...
                                      theta_min,eta,delta,gamma_theta,gamma_phi,s_phi,s_theta)

            x_k1 = [];
            % solve soc Q-SDP with adapted constraint
            [x_star_soc,skip_soc] = solve_Q_SDP_soc(obj,x_k,x_star,p0,Bk,args);
            
            % skip_soc means the Q-sdp might not produce a feasible
            % solution
            if ~skip_soc
                        
                    dk     = x_star     - x_k;
                    dk_soc = x_star_soc - x_k;
                    
                    % use corected search direction and check if adapted point is acceptable to filter
                    dk_corr = dk + dk_soc;

                     % check filter acceptance
                    [x_k1_soc,theta_x_k1_soc,f_x_k1_soc,filter_Acceptance] = check_filter_acceptance(obj,filter,alpha,x_k,dk_corr,p0,args);


                    if filter_Acceptance % acceptable to filter

                        % check sufficient decrease
                        [suffDecrease_flag,~,~,filter] = chechSuffDecrease(obj,alpha,x_star_soc,x_k,p0,theta_xk,theta_x_k1_soc, f_x_k1_soc,f_xk, ...
                                                                             theta_min,eta,delta,gamma_theta,gamma_phi,s_phi,s_theta,filter);
    
    
                        if suffDecrease_flag
                            % soc solution becomes solution
                            x_k1 = x_k1_soc;
                            
                        end

                    end
      
            end
end % end of function