function min_dtry = bisection_minimization(obj, xk, xk1, sol)
    % Bisection algorithm to minimize the given function
    
    % Define the function to minimize
    f = @(dtry) double(obj.Merit(xk.*(1-dtry) + xk1.*dtry, sol{5}));
    
    % Set the tolerance
    tolerance = 1e-6;
    
    % Initialize bounds (linesearch 0 <= d <= 1)
    lower_bound = 0;
    upper_bound = 1;
    
    % Start bisection loop
    while (upper_bound - lower_bound) > tolerance
        % Calculate mid-point
        mid_point = (lower_bound + upper_bound) / 2;
        
        % Calculate function values at end points and mid-point
        f_lower = f(lower_bound);
        f_mid   = f(mid_point);
        f_upper = f(upper_bound);
        
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
