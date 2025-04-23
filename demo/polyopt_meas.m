% Moment-SOS alternative

% indeterminate variable
x = casos.Indeterminates('x');
% some polynomial
f = x^4 + 10*x;

% create measure \mu with real sequence of moments y=(y_\alpha)
% interpret as a univariate polynomial 
deg = f.maxdeg;
s = casos.PS.sym('q',monomials(x,0:deg));

% get the cost function
cost = s.dot(f); % inner product between a polynomial and a measure

% get the decision variables and set the Hankel matrix ( for the
% multivariate case it is the moment matrix)
dec = [s.list_of_coeffs{:}];
dec = dec(:);

m = 0.5*deg+1;
H = casadi.SX(m, m);  % Initialize an empty symbolic matrix

% Fill the Hankel matrix using a loop
for i = 1:m
    for j = 1:m
        H(i, j) = dec(i + j - 1); % Hankel indexing rule
    end
end

% cone definition
opts = struct('Kc', struct('psd',m), 'Kx', struct('lin', length(dec)));

% build casos problem
meas = struct('x', dec, 'f', cost, 'g', H(:));

% solve by relaxation to SDP
D = casos.sdpsol('S','mosek',meas,opts);

% evaluate
warning('off', 'all');  % Turn off all warnings
sol = D('lbx', [1; -inf(length(dec)-1,1)] , 'ubx', [1; inf(length(dec)-1,1)] );

fprintf('Minimum is %g.\n', full(sol.f))
