function [xi_k1, dual_k1, Filter, Linesearch_stall,alpha,predictedcost,predictedConVio]   = lineSearchFilterSimple(obj,p0,Filter,xi_k,xi_plus,dual_k,dual_plus)   

%% backtracking line search filter (simple)
        alpha            = 1;
        rho              = 0.8;
        FilterAccept     = 0;  
        maxIter          = 0;
        Linesearch_stall = 0;
		
        while ~FilterAccept 
                
            xi_k1   = xi_k   + alpha*(xi_plus-xi_k);
            dual_k1 = dual_k + alpha*(dual_plus-dual_k);

            % Check predicted cost and predicted constraint violation
            solPara_proj = obj.projConPara('p',[xi_k1;p0]);

            predictedConVio = (double(solPara_proj.f));

            predictedcost   = double(obj.cost_fun(xi_k1));
        
            % setup logical array; 
            dom_logi_arr = Filter <= [predictedcost, predictedConVio]; 
              
            % if both entries are 1 than domination list entry dominates
            % current iterate
            if sum(dom_logi_arr,2) < 2

                % Remove all pairs in filter that are dominated by the current pair; 
                % if both entries are zero means both are smaller than list entry
                Filter(sum(dom_logi_arr,2) == 0 ,:) = [];
                
                % add new accepted pair from current iteration
                Filter       = vertcat(Filter, [predictedcost, predictedConVio]);

                FilterAccept = true;

            else
                alpha = alpha*rho;
                FilterAccept = false;
            end
            

            maxIter = maxIter  + 1;
            
            if maxIter >= 20 || alpha <= 0.001
                
                Linesearch_stall = 1;
                break
            end
        end

end
