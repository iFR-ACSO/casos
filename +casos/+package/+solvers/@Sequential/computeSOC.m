function [xi_k1,alpha_k] = computeSOC(obj,args_solver,p0,xi_k,xi_plus,H,s_phi,s_theta,eta,gamma_phi,gamma_theta,Filter,old_conVio,tau,alpha_k)
                    
                    % To-Do add constraint violation decrease by factor
                
                    pmax            = 4;
                    p               = 1;
                    alpha_soc       = alpha_k; 

                     gamma_theta = 0.01;
                     gamma_phi  =  1;

                    % initializes
                    d_corr_full = xi_plus- xi_k;
                    % orig
                    % SOC
                    while p < pmax % perform for maximum number of SOC iterations
                        
                        % compute Q-SDP subproblem with SOC correction term
                        d_corr = solveSOC_subSDP(obj,args_solver,p0,xi_k,d_corr_full,H); 
                        
                        % compute new search direction
                        d_corr_full = d_corr_full + d_corr;
                        
                        % compute new trial point based on corrected
                        % direction
                        x_soc_trial = xi_k + alpha_soc*d_corr_full;

                        % check if correct trial point is acceptable to
                        % filter
                        [filter_accept_soc,Filter,predictedConVio,predictedcost]  = checkFilterAcceptance(obj,Filter,x_soc_trial,p0,gamma_theta,gamma_phi);
                          
                       % if acceptable check if sufficient decrease conditions are fulfilled
                        if filter_accept_soc
                                
                            accept_soc = compSuffDecreCond_SOC(obj,...
                                                           p0,...
                                                           x_soc_trial,...
                                                           xi_plus,...
                                                           xi_k,...
                                                           alpha_k,...  % not sure about this; is alpha_k from outer-linesearch
                                                           alpha_soc,...
                                                           old_conVio,...
                                                           predictedConVio,...
                                                           predictedcost,...
                                                           s_phi,...
                                                           s_theta,...
                                                           eta,...
                                                           gamma_phi,...
                                                           gamma_theta,...
                                                           Filter);


                            if accept_soc 
              
                                % if fulfilled set alpha SOC to alpha
                                % lineserach and accept step
                                % alpha_k = alpha_soc;
                                xi_k1   = x_soc_trial;
                                return
                            end

                        else
                            % not acceptable to filter; adjust alpha from
                            % lineserach
                            alpha_k = alpha_k*tau;
                            % reduceAlphaLinesearch = 1;
                            xi_k1 = [];
                            return
                             
                            % alpha_soc = alpha_soc*tau ;
                        end
          

                        % counter
                        p = p + 1;

                        old_conVio = predictedConVio;

                    end % end SOC-while loop

                    if p == pmax

                        alpha_k = alpha_k*tau;
                        xi_k1 = [];
                    end

end