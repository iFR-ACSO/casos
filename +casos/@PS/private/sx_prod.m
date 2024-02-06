function b = sx_prod(a,dim)
% Product of SX matrix elements.

    if isempty(a)
        % product of empty matrix
        b = casadi.SX(size(a));
        return

    elseif size(a,dim) == 1
        % nothing to do
        b = a;
        return

    elseif dim < 2
        % fall back to second dimension
        b = sx_prod(a',2)'; 
        return
    end
    
    % else
    [m,n] = size(a);

    % create function f: (X,U) -> X.*U
    X = casadi.SX.sym('X',m,1);
    U = casadi.SX.sym('U',m,1);
    f = casadi.Function('f',{X U},{X.*U});
    
    % fold to F: (X0,[U1...Un]) -> f(f(...f(X0,U1),...),Un)
    F = fold(f,n);
    % evaluate matrix product
    b = F(1,a);
end
