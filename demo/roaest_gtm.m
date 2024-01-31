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
opts.Kx = struct('s', 1);
opts.Kc = struct('s', 1);
opts.error_on_fail = 0;
S1 = casos.qcsossol('S1','bisection',sos1,opts);

% solver 2: beta-step
sos2 = struct('x',s2,'f',-b,'p',[V;g]);
sos2.('g') = s2*(p-b)+g-V;

% states + constraint are SOS cones
opts.Kx = struct('s', 1);
opts.Kc = struct('s', 1);

S2 = casos.qcsossol('S2','bisection',sos2,opts);

% solver 3: V-step
s1_sym = casos.PS.sym('s1',basis(s1));
s2_sym = casos.PS.sym('s2',basis(s2));

% s1 = casos.PS.sym()
Vlb = casos.PS(basis(V),-inf);
Vub = casos.PS(basis(V),+inf);

sos3 = struct('x',V,'p',[b,g,s1_sym,s2_sym]);
sos3.('g') = [V-l; s2_sym*(p-b)+g-V; s1_sym*(V-g)-nabla(V,x)*f-l];

opts = struct;
opts.Kx = struct('s', 0, 'l', 1); 
opts.Kc = struct('s', 3);

S3 = casos.sossol('S','sedumi',sos3,opts);

%% gamma-beta-V-iteration

for iter = 1:3

    % gamma step
    sol1 = S1('p',Vval);

    gval = double(-sol1.f);
    s1val = sol1.x;

    % beta step
    sol2 = S2('p',[Vval;gval]);

    bval = -sol2.f;
    s2val = sol2.x;

    % V-step
    sol3 = S3('p',[bval,gval,s1val,s2val],'lbx',Vlb,'ubx',Vub);


    Vval = sol3.x;

end


d = ([
          convvel(20, 'm/s', 'm/s')  %range.tas.lebesgue.get('ft/s')
          convang(20, 'deg', 'rad')  %range.gamma.lebesgue.get('rad')
          convang(50, 'deg', 'rad')  %range.qhat.lebesgue.get('rad')
          convang(20, 'deg', 'rad')  %range.alpha.lebesgue.get('rad')
]);
% d = ones(4,1);

D = diag(d)^-1;

V = subs(Vval,x,D*x);
V1 = subs(V,[x(1);x(4)],[0 0]');

p = subs(p,x,D*x);
p1 = subs(p,[x(1);x(4)],[0 0]');

figure(1)
pcontour(V1, double(gval), [-1 1 -4 4], 'b-');
hold on
pcontour(p1, double(bval), [-1 1 -4 4], 'r--');
