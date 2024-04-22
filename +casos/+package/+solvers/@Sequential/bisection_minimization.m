function min_dtry = bisection_minimization(obj,p0, xk, xk1, dual_plus )
    % Bisection algorithm to minimize the given function
    
    % Define the function to minimize
    % f = @(dtry) double(obj.Merit(xk.*(1-dtry) + xk1.*dtry, dual_plus,p0));
    
    % Set the tolerance
    tolerance = 1e-4;
    
    % Initialize bounds (linesearch 0 <= d <= 1)
    lower_bound = 0.1;
    upper_bound = 1;
    
    % Start bisection loop
    while (upper_bound - lower_bound) > tolerance
        % Calculate mid-point
        mid_point = (lower_bound + upper_bound) / 2;
        
        % Calculate function values at end points and mid-point
        f_lower = double(obj.Merit(xk.*(1-lower_bound) + xk1.*lower_bound, dual_plus,p0));% f(lower_bound);
        f_mid   = double(obj.Merit(xk.*(1-mid_point) + xk1.*mid_point, dual_plus,p0));% f(mid_point);
        f_upper =double(obj.Merit(xk.*(1-upper_bound) + xk1.*upper_bound, dual_plus,p0));%  f(upper_bound);
        
        % Check signs of function values
        if (f_mid * f_lower) < 0
            % Change upper bound if signs are different
            upper_bound = mid_point;
        elseif (f_mid * f_upper) < 0
            % Change lower bound if signs are different
            lower_bound = mid_point;
        else
            % If signs are the same, reduce the interval
            lower_bound = mid_point;
            upper_bound = mid_point;
        end
    end
    
    % Return the midpoint as the minimum dtry
    min_dtry = (lower_bound + upper_bound) / 2;
end
