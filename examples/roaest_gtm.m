% Estimate region of the Generic Transport Model
% See Chakraborty et al. 2011 (CEP) for details.

% system states
x = casos.PS('x',4,1);

% Polynomial Dynamics
load GTM_scaled_dyn.mat

f = gtmdyn(x(1),x(2),x(3),x(4));

% shape function
p = x'*x*1e2;

% initial Lyapunov functionn
P = [395.382671347059	-23.0032507976836	3.16965275615691	29.2992065909380
    -23.0032507976836	90.4764915483638	-16.1191428789579	-132.594376986429
    3.16965275615691	-16.1191428789579	3.44002648214490	24.2292666551058
    29.2992065909380	-132.594376986429	24.2292666551058	202.114797577027];
Vval = x'*P*x;

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

% options
opts = struct('sossol','sedumi');

%% setup solver

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

%% gamma-beta-V-iteration

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

end
