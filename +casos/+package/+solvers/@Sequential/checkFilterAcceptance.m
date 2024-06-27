function [filter_accept, Filter,predictedConVio,predictedcost] = checkFilterAcceptance(obj,Filter,x_trial,p0,gamma_theta,gamma_phi)

  % % Check predicted cost and predicted constraint violation
  solPara_proj = obj.projConPara('p',[x_trial;p0]);

  predictedConVio = double(solPara_proj.f);


               %    conVio_fun = to_function(obj.constraintFun(x_trial,p0));
               % 
               % conVioVal = full(conVio_fun(obj.samplingPoint_proj));
               %  if min(min(conVioVal)) > 0
               %      predictedConVio = 0;
               %  else
               %      predictedConVio = abs(min(min(conVioVal)));
               %  end
            
  predictedcost   = double(obj.cost_fun(x_trial));
  
  % From Fletcher
  % beta   = 0.9;

  Filter_adapt = [Filter(:,1) - gamma_phi*Filter(:,2), Filter(:,2)*(1-gamma_theta)];
    
  % For acceptance one cost or violation must be smaller
  % Example comparison one list entry with current trial point:
  % Output: 1) 0 0 --> Filter dominates i.e. entry is smaller
  %         2) 0 1 --> cost larger, but violation smaller 
  %         3) 1 0 --> cost smaller, but violation larger
  %         4) 1 1 --> both are smaller
  % --> case 2) to case 4): accept to filter and augment
  dom_logi_arr = Filter_adapt > [predictedcost, predictedConVio]; 
                          
  % if both entries are 1 than  list entry dominates current iterate
  if sum(dom_logi_arr,2) ~= 0
      filter_accept = 1;
  else
      filter_accept = 0;
  end

end % end of function