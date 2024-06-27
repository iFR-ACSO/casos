function d_corr = solveSOC_subSDP(obj,args_solver,p0,xi_k,xi_plus,H)

                        args_solver{2} = [p0; xi_k;xi_plus;H(:)];

                         % evaluate subproblem of SOC
                        sol_soc = call(obj.solverSOC, args_solver);
                        

                        switch (obj.solverSOC.stats.UNIFIED_RETURN_STATUS)
                            case 'SOLVER_RET_SUCCESS'
                                
               
                                d_corr = sol_soc{1};
                                
                
                            case 'SOLVER_RET_UNKNOWN'
                                 % problem seems feasible, but but solution is not
                
                                 d_corr = sol_soc{1};
                                
                                
                            case 'SOLVER_RET_INFEASIBLE'
                                
                                 error('SOC infeasible')
                
                            otherwise   % problem status unknown or ill-posed
                                
                               d_corr = sol_soc{1};
                            
                        end


end