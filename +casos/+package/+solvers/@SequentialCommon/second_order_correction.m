%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Short Description:
%
% This function implements a simple second-order correction step. An
% adapted quadratic subproblem is solved to correct the search direction.
% For the new search direction filter acceptance and sufficient decrease is
% checked. If not succesful, than we return to original probelm
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_k1,suffDecrease_flag,f_type,amijo] = second_order_correction(obj, ...
    x_k, ...
    x_star, ...
    p0, ...
    Bk, ...
    args, ...
    filter, ...
    alpha, ...
    theta_xk, ...
    f_xk, ...
    dkl, ...
    L_k1, ...
    L_k, ...
    dual_k)


% solve soc Q-SDP with adapted constraint
[x_star_soc,dual_star_soc,skip_soc] = solve_Q_SDP_soc(obj,x_k,x_star,p0,Bk,args);

% skip_soc means the Q-sdp might not produce a feasible
% solution
if ~skip_soc

    % actual search direction
    dk     = x_star     - x_k;
    % new search direction
    dk_soc = x_star_soc - x_k;

    % use corected search direction and check if adapted point is acceptable to filter
    dk_corr = dk + dk_soc;
    % correct dual direction
    dkl_corr = dkl + (dual_star_soc-dual_k);

    % check filter acceptance for corrected search direction
    [x_k1_soc,theta_x_k1_soc,f_x_k1_soc,filter_Acceptance] = check_filter_acceptance(obj,filter,alpha,x_k,dual_k,dk_corr,dkl_corr,p0,args);

    if filter_Acceptance
        % check sufficient decrease for corrected search direction
        [suffDecrease_flag,f_type,amijo] = check_suffDecrease(obj,alpha,x_star_soc,x_k,p0,theta_xk,theta_x_k1_soc, f_x_k1_soc,f_xk,L_k1,L_k,dual_k,filter);

        if suffDecrease_flag
            % soc solution becomes solution of the actual problem
            x_k1 = x_k1_soc;
        else
            x_k1 = x_k;
            suffDecrease_flag = 0;
            f_type = 0;
            amijo = 0;

        end
    else
        x_k1 = x_k;
        suffDecrease_flag = 0;
        f_type = 0;
        amijo = 0;
    end
else
    x_k1 = x_k;
    suffDecrease_flag = 0;
    f_type = 0;
    amijo = 0;

end

end % end of function