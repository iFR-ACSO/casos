% Estimate region of attraction by V-s-iteration and bisection.

% system states
x = casos.Indeterminates('x',2);

% system dynamics
f = [-x(2); x(1) + (x(1)^2 - 1)*x(2)];

% Lyapunov function candidate
Vval = 1.5*x(1)^2 - x(1)*x(2) + x(2)^2;
p = x'*x;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2));

% SOS multiplier
s1 = casos.PS.sym('s1',monomials(x,1),'gram');
s2 = casos.PS.sym('s2',monomials(x,0),'gram');

% enforce positivity
l = 1e-6*(x'*x);

% level of stability
g = casos.PS.sym('g');
b = casos.PS.sym('b');

% options
opts = struct('sossol','sedumi');

%% Setup solver
% solver 1: gamma-step
sos1 = struct('x',s1,'f',-g,'p',V);
sos1.('g') = s1*(V-g)-nabla(V,x)*f-l;

% states + constraint are SOS cones
opts.Kx = struct('sos', 1);
opts.Kc = struct('sos', 1);

S1 = casos.qcsossol('S1','bisection',sos1,opts);

% solver 2: beta-step
sos2 = struct('x',s2,'f',-b,'p',[V;g]);
sos2.('g') = s2*(p-b)+g-V;

% states + constraint are SOS cones
opts.Kx = struct('sos', 1);
opts.Kc = struct('sos', 1);

S2 = casos.qcsossol('S2','bisection',sos2,opts);

% solver 3: V-step
sos3 = struct('x',V,'p',[b,g,s1,s2]);
sos3.('g') = [V-l; s2*(p-b)+g-V; s1*(V-g)-nabla(V,x)*f-l];

opts = struct;
opts.Kx = struct('sos', 0, 'lin', 1); 
opts.Kc = struct('sos', 3);

S3 = casos.sossol('S','sedumi',sos3,opts);

%% V-s-iteration
for iter = 1:10

    % gamma step
    sol1 = S1('p',Vval);

    gval = -sol1.f;
    s1val = sol1.x;

    % beta step
    sol2 = S2('p',[Vval;gval]);

    bval = -sol2.f;
    s2val = sol2.x;

    % V-step
    sol3 = S3('p',[bval,gval,s1val,s2val]);

    Vval = sol3.x;

    fprintf('Iteration %d: b = %g, g = %g.\n',iter,full(bval),full(gval));
end
