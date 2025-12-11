function [map, Qvar_G, Qcon_G] = process_reorder(Qvar_G, Qcon_G, ...
                                                 Ksdp_x_s, Ksdp_g_s,  ...
                                                 Qvar_l, Ns, Ms, Nds, Mds)
        % Reorder variables to group by cone type in SDP formulation.
        %
        % This function creates a permutation matrix to reorder SDP variables so that
        % variables are grouped by cone type: linear >> PSD >> DD >> SDD.
        %
        % Inputs:
        %   Qvar_G    - Gram matrix variables for SOS variables
        %   Qcon_G    - Gram matrix variables for SOS constraints  
        %   Ksdp_x_s  - Matrix dimensions for variable cones
        %   Ksdp_g_s  - Matrix dimensions for constraint cones
        %   Qvar_l    - Linear variables
        %   Ns        - Number of SOS variable cones
        %   Ms        - Number of SOS constraint cones
        %   Nds       - Number of DSOS variable cones
        %   Mds       - Number of DSOS constraint cones
        %
        % Outputs:
        %   map       - Permutation matrix for reordering
        %   Qvar_G    - Reordered variable Gram matrices
        %   Qcon_G    - Reordered constraint Gram matrices

        % PSD variable cumulative sizes
        cum_vars = [0; cumsum(Ksdp_x_s.^2)];  

        % PSD constraint cumulative sizes
        cum_cons = [0; cumsum(Ksdp_g_s.^2)];  
        
        % linear variables offset
        base_lin = numel(Qvar_l);    
        % variable group offset
        base_var = numel(Qvar_G);   
        % constraint group offset
        base_con = numel(Qcon_G);         
    
        % linear variable indices
        linear_indices = 1:base_lin;
        
        % PSD sections
        psd_var_indices = arrayfun(@(i) base_lin + (cum_vars(i)+1:cum_vars(i+1)), 1:Ns, 'UniformOutput', false);
        psd_con_indices = arrayfun(@(i) base_lin + base_var + (cum_cons(i)+1:cum_cons(i+1)), 1:Ms, 'UniformOutput', false);
        
        % DD sections
        dd_var_indices = arrayfun(@(i) base_lin + (cum_vars(Ns+i)+1:cum_vars(Ns+i+1)), 1:Nds, 'UniformOutput', false);
        dd_con_indices = arrayfun(@(i) base_lin + base_var + (cum_cons(Ms+i)+1:cum_cons(Ms+i+1)), 1:Mds, 'UniformOutput', false);
        
        % SDD sections
        sdd_var_indices = base_lin + (cum_vars(Ns+Nds+1)+1:cum_vars(end));
        sdd_con_indices = base_lin + base_var + (cum_cons(Ms+Mds+1)+1:cum_cons(end));
        
        % combine all indices
        all_indices = [ ...
            linear_indices, ...
            [psd_var_indices{:}], [psd_con_indices{:}], ...
            [dd_var_indices{:}], [dd_con_indices{:}], ...
            sdd_var_indices, sdd_con_indices ...
        ];
        
        % build permutation matrix
        total_size = base_lin + base_var + base_con;
        if numel(all_indices) ~= total_size
            warning('Index length mismatch: expected %d, got %d', total_size, numel(all_indices));
        end
        
        map = sparse(1:total_size, all_indices, 1, total_size, total_size);

        % new Qvar_sdp and Qcon_G
        sdp_x = map*[Qvar_l; Qvar_G; Qcon_G]; % reorder the sdp

        Qvar_G = sdp_x([linear_indices, psd_var_indices{:}, dd_var_indices{:}, sdd_var_indices]);
        Qcon_G = sdp_x([psd_con_indices{:}, dd_con_indices{:}, sdd_con_indices]);
        
    end