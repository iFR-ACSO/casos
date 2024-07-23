classdef Filter 
    

    properties (Access=private)
        filter = [];
        opts

    end

    properties (Constant,Access=protected)
     filter_options = [
        {'maxConVio_0', 'Maximum allowed constraint violation.'
         'maxConVio_0_multiplier', 'Multiplier for filter initialization.'
         'gamma_theta', 'Filter envelope parameter for constraint violation.'
         'gamma_phi', 'Filter envelope parameter for cost.'
         'eta', 'Scaling parameter for Amijo-condition'
         'minConVio', 'Minimum threshold in constraint violation'
         's_theta' , 'Scaling parameter for sufficient decrease in constraint violation' 
         's_phi' , 'Scaling parameter for sufficient decrease in cost' 
         }
    ];

    end

    properties (SetAccess=private)
        class_name = 'Filter';
    end


    methods
        % constructor
        function obj = Filter(specOps)
            
             % default options
            if ~isfield(specOps,'maxConVio_0'),             obj.opts.maxConVio_0            = 10000000000;  end
            if ~isfield(specOps,'maxConVio_0_multiplier'),  obj.opts.maxConVio_0_multiplier = 100;          end
            if ~isfield(specOps,'gamma_theta'),             obj.opts.gamma_theta            = 0.01;         end
            if ~isfield(specOps,'gamma_phi'),               obj.opts.gamma_phi              = 1;            end
            if ~isfield(specOps,'eta'),                     obj.opts.eta                    = 1e-4;         end
            if ~isfield(specOps,'s_phi'),                   obj.opts.s_phi                  = 2.3;          end
            if ~isfield(specOps,'s_theta'),                 obj.opts.s_theta                = 1.1;          end
            if ~isfield(specOps,'minConVio'),               obj.opts.minConVio              = 1e-4;         end
  
        end
        
        % method to initialize filter before the first iteration
        function obj = initializeFilter(obj,conVio_xi0,cost0)
            

             % compute most right entry of filter
             conVio_0   = max(obj.opts.maxConVio_0, obj.opts.maxConVio_0_multiplier*conVio_xi0);

             % compute most right entry of filter
             % conVio_0   = obj.opts.maxConVio_0_multiplier*conVio_xi0;
            
             % minimum constraint violation for sufficient decrease
             obj.opts.minConVio =  1e-3;
                
             % initialize filter; basically an array
             obj.filter             = [cost0,conVio_0];

        end
        
        % function called during online execution
        function [obj,goto_SOC,accept_new_iter] = updateFilter(obj,...
                                                               alpha_k,...
                                                               searchDir, ...
                                                               curr_cost, ....
                                                               curr_conVio,....
                                                               new_cost,...
                                                               new_conVio,....
                                                               nabla_x_cost)

            
            % check if current iterate is acceptable to filter
            [obj,FilterAcceptFlag] = checkFilterAcceptance(obj,new_cost,new_conVio);

            if FilterAcceptFlag
                    

                % check sufficient decrease
                [obj,suffDecreaseFlag] = checkSuffDecrease(obj,...
                                                            alpha_k,...
                                                            searchDir, ...
                                                            curr_cost, ....
                                                            curr_conVio,....
                                                            new_cost,...
                                                            new_conVio,....
                                                            nabla_x_cost);
                
                
                if suffDecreaseFlag
                    % accept new iterate; leave linesearch
                    goto_SOC            = 0;
                    accept_new_iter     = 1; 
                else
                   % go to second-order-correction
                   accept_new_iter      = 0;
                   goto_SOC             = 1;
                end


            else
                % decrease step-length
                accept_new_iter     = 0;
                goto_SOC            = 0;
            end

        end % --- end of method updateFilter ---


    end

    methods (Access =private)
        % check if a new iterate can be accepted to filter or not   
        function [obj,FilterAcceptFlag] = checkFilterAcceptance(obj,new_cost,new_conVio)
                
              % add envelope
              Filter_adapt = [obj.filter(:,1) - obj.opts.gamma_phi*obj.filter(:,2), ... 
                              obj.filter(:,2)*(1-obj.opts.gamma_theta)];
    
              % For acceptance one cost or violation must be smaller
              % Example comparison one list entry with current trial point:
              % Output: 1) 0 0 --> Filter dominates i.e. entries are smaller
              %         2) 0 1 --> cost larger, but violation smaller 
              %         3) 1 0 --> cost smaller, but violation larger
              %         4) 1 1 --> both are smaller
              % --> case 2) to case 4): accept to filter 
              dom_logi_arr = Filter_adapt > [new_cost,new_conVio]; 
                                      
              % if both entries are 1 than  list entry dominates current iterate
              if sum(dom_logi_arr,2) ~= 0  % row-wise sum must be unequal to 0
                  FilterAcceptFlag = 1;
                     % augment filter
                   obj.filter = vertcat(obj.filter, [new_cost,new_conVio]);
              else
                  FilterAcceptFlag = 0;
              end


        end % --- end of function checkFilterAcceptance ---

        % check for sufficient decrease
        function [obj,suffDecreaseFlag] = checkSuffDecrease(obj,...
                                                            alpha_k,...
                                                            searchDir, ...
                                                            curr_cost, ....
                                                            curr_conVio,....
                                                            new_cost,...
                                                            new_conVio,....
                                                            nabla_x_cost)

            
            % minimum decrease in constraint violation
            suffDecrease_conVio = curr_conVio < obj.opts.minConVio;
            
            % sufficient decrease in cost; descent direction and decrease
            suffDecrease_cost   = nabla_x_cost'*searchDir < 0 &&  ...
               alpha_k*( - nabla_x_cost'*searchDir )^obj.opts.s_phi > curr_conVio^obj.opts.s_theta;
            

                % check for f-type switching
                if suffDecrease_conVio && suffDecrease_cost
                        
                        % check if Amijo-condition is fulfilled
                        if new_cost <= curr_cost + alpha_k*obj.opts.eta*nabla_x_cost'*searchDir
                            suffDecreaseFlag = 1;
                        else
                            suffDecreaseFlag = 0;
                        end

                else

                   % augment filter
                   % obj.filter = vertcat(obj.filter, [curr_cost, curr_conVio]);
                    
                   % check if the cost or constraint violation of new iterate is at least as good as the
                   % current iterate
                   if new_conVio    <= (1-obj.opts.gamma_theta)*curr_conVio ||...
                      new_cost      <= curr_cost - obj.opts.gamma_phi*curr_conVio

                      suffDecreaseFlag = 1;

                   else

                       suffDecreaseFlag = 0;

                   end

                end

        end % --- end of function checkSuffDecrease ---

    end % --- end of private methods ---

end % --- end of class filter ---