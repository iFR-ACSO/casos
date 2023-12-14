% Validate stability with LF candidate.
clear
clc

% system states
x = casos.PS('x',2,1);

% system dynamics
f = [-x(2); x(1) + (x(1)^2 - 1)*x(2)];

% Lyapunov function candidate
Vval = 1.5*x(1)^2 - x(1)*x(2) + x(2)^2;
p = Vval*5;

% derivative w.r.t. time
Vvaldot = nabla(Vval,x)*f;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,1:4));

% SOS multiplier
s1 = casos.PS.sym('s1',monomials(x,1:4),'gram');
s2 = casos.PS.sym('s2',monomials(x,1:4),'gram');

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
opts.error_on_fail = false;

%% setup solver

% solver 1: gamma-step
sos1 = struct('x',s1,'f',-g,'g',s1*(V-g)-nabla(V,x)*f-l,'p',V);

% states + constraint are SOS cones
opts.Kx.s = 1; 
opts.Kc.s = 1;

% S1 = casos.sossol('S','sedumi',sos1,opts);
S1 = casos.qcsossol('S1','bisection',sos1,opts);

% solver 2: beta-step
sos2 = struct('x',s2,'f',-b,'g', s2*(p-b)+g-V,'p',[V;g]);

% states + constraint are SOS cones
opts.Kx.s = 1; 
opts.Kc.s = 1;

% S2 = casos.sossol('S','sedumi',sos2,opts);
S2 = casos.qcsossol('S2','bisection',sos2,opts);

% solver 3: V-step
[~,R,~] =poly2basis(s1);
s1_sym = casos.PS.sym('s1',R);
[~,R,~] =poly2basis(s2);
s2_sym = casos.PS.sym('s2',R);

% s1 = casos.PS.sym()
Vlb = casos.PS(basis(V),-inf);
Vub = casos.PS(basis(V),+inf);

sos3 = struct('x',V,...
              'g',[V-l; s2_sym*(p-b)+g-V;s1_sym*(V-g)-nabla(V,x)*f-l],...
              'p',[b,g,s1_sym,s2_sym]);

opts.Kx.s = 0; 
opts.Kx.l = 1; 
opts.Kc.s = 3;

S3 = casos.sossol('S','sedumi',sos3,opts);

%% gamma-beta-V-iteration

tic 

% profile on

for iter = 1:10


    sol1 = S1('p',Vval);

    gval = double(-sol1.f);
    s1val = sol1.x;

    sol2 = S2('p',[Vval;gval]);

   bval = -sol2.f;
   s2val = sol2.x;


 % V-step
 sol3 = S3('p',[bval,gval,s1val,s2val],'lbx',Vlb,'ubx',Vub);


    switch (S3.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS',... 
                Vval = sol3.x;

        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}
            disp('V step infeasible')
            break
        otherwise, error('Failed.')
    end
end

% profile viewer
toc



