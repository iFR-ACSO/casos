% Estimate region of attraction by V-s-iteration.

% system states
x = casos.PS('x',2,1);

% system dynamics
f = [-x(2); x(1) + (x(1)^2 - 1)*x(2)];

% Lyapunov function candidate
Vval = 1.5*x(1)^2 - x(1)*x(2) + x(2)^2;
p = Vval*5;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2));

% SOS multiplier
s1 = casos.PS.sym('s1',monomials(x,0:2),'gram');
s2 = casos.PS.sym('s2',monomials(x,0:2),'gram');

% enforce positivity
l = 1e-6*(x'*x);

% level of stability
g = casos.PS.sym('g');
b = casos.PS.sym('b');

% ignore infeasibility
opts = struct('error_on_fail', false);

%% setup solver

% solver 1: gamma-step
sos1 = struct('x',s1,'p',[g;V]);
sos1.('g') = s1*(V-g)-nabla(V,x)*f-l;

% states + constraint are SOS cones
opts.Kx = struct('s', 1);
opts.Kc = struct('s', 1);

S1 = casos.sossol('S','sedumi',sos1,opts);

% solver 2: beta-step
sos2 = struct('x',s2,'p',[g;b;V]);
sos2.('g') = s2*(p-b)+g-V;

% states + constraint are SOS cones
opts.Kx = struct('s', 1);
opts.Kc = struct('s', 1);

S2 = casos.sossol('S','sedumi',sos2,opts);

% solver 3: V-step
s1_sym = casos.PS.sym('s1',basis(s1));
s2_sym = casos.PS.sym('s2',basis(s2));

% s1 = casos.PS.sym()
Vlb = casos.PS(basis(V),-inf);
Vub = casos.PS(basis(V),+inf);

sos3 = struct('x',V,'p',[b,g,s1_sym,s2_sym]);
sos3.('g') = [V-l; s2_sym*(p-b)+g-V; s1_sym*(V-g)-nabla(V,x)*f-l];

opts.Kx = struct('s', 0, 'l', 1); 
opts.Kc = struct('s', 3);

S3 = casos.sossol('S','sedumi',sos3,opts);

%% V-s-iteration

for iter = 1:10

% gamma step: find largest possible stable level set
lb = 0; 
ub = 1000;

% bisection
while (ub - lb > 1e-1)
    s1val = [];
    s2val = [];
    gtry = (lb+ub)/2;

    % evaluate parametrized SOS problem
    sol1 = S1('p',[gtry;Vval]);

    switch (S1.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS'
            s1val = sol1.x;
            lb = gtry;
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}
            ub = gtry;
        otherwise, error('Failed.')
    end
end

if ~isempty(s1val)
    fprintf('Maximal stable level set is %g.\n', lb);
else
    disp('Not feasible in gamma step!')
    break
end
gval = lb;


% beta step: find largest possible stable level set
lb = 0; 
ub = 1000;

% bisection
while (ub - lb > 1e-1)
    btry = (lb+ub)/2;

    % evaluate parametrized SOS problem
    sol2 = S2('p',[gval;btry;Vval]);

    switch (S2.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS'
            lb = btry;
            s2val = sol2.x;
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}
            ub = btry;
        otherwise, error('Failed.')
    end
end
if ~isempty(s2val)
    fprintf('Maximal shape function level set is %g.\n', lb);
else
    disp('Not feasible in beta step!')
    break
end

bval = lb;


% V-step
sol3 = S3('p',[bval,gval,s1val,s2val],'lbx',Vlb,'ubx',Vub);

switch (S3.stats.UNIFIED_RETURN_STATUS)
    case 'SOLVER_RET_SUCCESS'
            Vval = sol3.x;

    case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}
        disp('V step infeasible')
        break
    otherwise, error('Failed.')
end

end
