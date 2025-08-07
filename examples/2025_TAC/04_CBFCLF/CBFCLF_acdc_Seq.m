%% ------------------------------------------------------------------------
%
%
%   Short Descirption:  Calculate an inner-estimate of the
%                       region-of-attraction for the longitudinal motion 
%                       of the Nasa Generic Transport Model. To increase
%                       the size of the sublevel set we try to minimize the
%                       squared distance to a defined set. Additionally, we
%                       synthesis a linear control law at the same time.
%                       State and control constraints are considered.
%
%   Reference: Modified problem from:
%              Chakraborty, Abhijit and Seiler, Peter and Balas, Gary J.,
%              Nonlinear region of attraction analysis for flight control 
%              verification and validation, Control Engineering Practice,
%              2011, doi: 10.1016/j.conengprac.2010.12.001
%           
%
%--------------------------------------------------------------------------

clc
clear
close all

%% system states
x = casos.Indeterminates('x',3,1);
u = casos.Indeterminates('u',2,1);


f0 = [ -0.05*x(1) - 57.9*x(2) + 0.00919*x(3);
       1710*x(1) + 314*x(3);
      -0.271*x(1) - 314*x(2) ];
            
G = [ 0.05 - 57.9*x(2),   -57.9*x(3);
      1710 + 1710*x(1),        0;
                 0,  1710 + 1710*x(1) ];



f = f0+G*u;

A = full(subs(nabla(f,x),[x;u],[zeros(3,1);zeros(2,1)]));
B = full(subs(nabla(f,u),[x;u],[zeros(3,1);zeros(2,1)]));


[K0,P0] = lqr(A,B,eye(3),eye(2));

% w1(x)
w1 = (1*x(1) + 0.3)^2 + (1*x(2)/20)^2 + (1*x(3)/20)^2 - 0.5^2;

% w2(x)
w2 = (1*x(1)/20)^2 + 1*x(2)^2 + 1*x(3)^2 - 1.2^2;

% operating region
r=(1*x(1)/0.8)^2+(1*x(2)/1.2)^2+(1*x(3)/1.2)^2-1.8;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials([x(1)^2 x(2)^2 x(3)^2 x(1)^4 x(2)^4 x(3)^4]));

% SOS multiplier
s1    = casos.PS.sym('s1',monomials(x,0));
s4    = casos.PS.sym('s4',monomials(x,0),[2 1]);
s3    = casos.PS.sym('s3',monomials(x,2),[3 1]);
B1    = casos.PS.sym('b1',monomials([x(1)^2 x(2)^2 x(3)^2 x(1)^4 x(2)^4 x(3)^4]));
B2    = casos.PS.sym('b2',monomials([x(1)^2 x(2)^2 x(3)^2 x(1)^4 x(2)^4 x(3)^4]));


p = casos.PS.sym('p',monomials(x,1),[2 1]); 
q = casos.PS.sym('q',monomials(x,0:2),[2 1]);

% enforce positivity
l = 1e-6*(x'*x);


%% setup solver
% options
opts = struct('sossol','mosek');
opts.verbose  = 1;
opts.max_iter = 100;

% cost
cost = dot(r-B1,r-B1) + dot(r-B2,r-B2) ;

sos = struct('x',[V;B1;B2;s1;s3;s4;p;q],... % decision variables
              'f',cost, ...                   % cost
              'p',[]);                        % parameter

% SOS constraints
sos.('g') = [
             s1;
             s4
             s3;
             V-l; 
             - nabla(V,x)*(s1*f0+G*p) - l + s3(1)*r;
             - nabla(B1,x)*(s1*f0+G*p)  - q(1)*B1 + s3(2)*r;
             B1 - s4(1)*w1;
             - nabla(B2,x)*(s1*f0+G*p)  - q(2)*B2 + s3(3)*r;
             B2 - s4(2)*w2
             ];

% states + constraint cones
opts.Kx      = struct('lin', length(sos.x));
opts.Kc      = struct('sos', length(sos.g));

% setup solver
S = casos.nlsossol('S1','sequential',sos,opts);

%% solve problem

% initial guess
V0  = x'*P0*x;
B10 = x'*x;
B20 = x'*x;
s10 = x'*x;
s30 = x'*x*ones(3,1);
s40 = x'*x*ones(2,1);
p0  = -K0*x;
q0  = x'*x*ones(2,1);

x0 = casos.PD([V0;B10;B20;s10;s30;s40;p0;q0]);

% solve problem
sol = S('x0' ,x0);

casos.postProcessSolver(S,true);

S.stats

S.stats.single_iterations{end}.conic
