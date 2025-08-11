function b = casadi_prod(a,dim)
% Product of Casadi matrix elements.

if isempty(a)
    % product of empty matrix
    b = zeros(size(a),'like',a);
    return

elseif size(a,dim) == 1
    % nothing to do
    b = a;
    return

elseif dim < 2
    % fall back to second dimension
    b = casadi_prod(a',2)'; 
    return
end

% else
[m,n] = size(a);

% create function f: (X,U) -> X.*U
X = casadi.MX.sym('X',m,1);
U = casadi.MX.sym('U',m,1);
f = casadi.Function('f',{X U},{X.*U});

% fold to F: (X0,[U1...Un]) -> f(f(...f(X0,U1),...),Un)
F = fold(f,n);
% evaluate matrix product
b = F(1,a);

end
