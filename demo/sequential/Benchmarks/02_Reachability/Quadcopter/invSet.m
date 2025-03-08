%--------------------------------------------------------------------------
% 
% Implementation of custom V-s-iteration for the GTM 4D reachability 
% problem in CaSoS. The quasi-convex solver is used to perform the 
% bisection.
%
% Model and implementation of constraints based on 
% https://github.com/heyinUCB/Backward-Reachability-Analysis-and-Control-Synthesis
% 
%--------------------------------------------------------------------------

close all
clear
clc


% system states 
x = casos.Indeterminates('x',4);
t = casos.Indeterminates('t');

x1 = x(1); % V 
x2 = x(2); % alpha
x3 = x(3); % q
x4 = x(4); % theta

% Polynomial Dynamics
d2r = pi/180;
% scale matrix
Dmax = diag([20 20*d2r 50*d2r 20*d2r]);

f2 = - 10.6560*x1^3 + 11.5309*x1^2*x3 + 7.8850*x1*x3^2 + 0.7972*x2^2*x3 ...
       + 0.8408*x2*x3*x4 + 21.0492*x3^3 + 0.4204*x3*x4^2 + 66.5225*x1 - 24.5110*x3;
f4 = 10.9955*x1^3 - 48.9151*x1^2*x3 - 6.4044*x1*x3^2 - 2.3955*x2^2*x3 ...
    - 1.5943*x2*x3*x4 - 51.9088*x3^3 - 0.7971*x3*x4^2 - 68.6419*x1 + 103.9783*x3;
f = [x2; f2; x4; f4];

g2 = -10.0959*x3^2 + 44.2521;
g4 =  37.8015*x3^2 - 83.9120;

g = [0; g2; 0; g4];

% terminal set
rT = x'*blkdiag(1/0.1^2, 1/0.35^2, 1/0.1^2, 1/0.35^2)*x - 1;
% r0 = x'*blkdiag(1/0.2^2, 1/0.7^2, 1/0.2^2, 1/0.7^2)*x - 4;


% initial storage function
Vval = 182.2102*x1^2 + 67.8981*x1*x2 + 314.3265*x1*x3 + 37.2705*x1*x4 + ...
    6.4123*x2^2 + 59.0528*x2*x3 + 7.0314*x2*x4 + 138.9343*x3^2 + ...
    32.5944*x3*x4 + 1.9521*x4^2;

r0 = Vval-1;

% control constraint
uM =  1;
um = -1;


% Lyapunov function candidate
V    = casos.PS.sym('v',monomials([x],0:4));
K    = casos.PS.sym('k',monomials([x],1));

% SOS multiplier
s2 = casos.PS.sym('s2',monomials([x],0:2));
s3 = casos.PS.sym('s3',monomials([x],0:4));
s4 = casos.PS.sym('s4',monomials(x,0:2));
s5 = casos.PS.sym('s5',monomials([x],0:2));
s6 = casos.PS.sym('s6',monomials([x],0:2));
s7 = casos.PS.sym('s7',monomials([x],0:2));
s8 = casos.PS.sym('s8',monomials([x],0:2));
s9 = casos.PS.sym('s9',monomials(x,0:2));

b    = casos.PS.sym('b');

% options
opts.verbose = 1;
opts.max_iter = 150;


b = 1;
%% setup solver
sos1 = struct('x',[V;K;s3;s4;s5;s6;s7], ... % dec.var
              'f',dot( r0-(V-b), r0-(V-b) ), ... % cost function for bisection
              'p',[]);     % parameter

% constraint
sos1.('g') = [
             s4;
             s5;
             s7;
             s3*(V - b) -  nabla(V, x)*(f + g*K) ; ...
             s4*(V - b) - rT;...
             uM - K + s5*(V - b); ...
             K - um + s7*(V - b )
             ];


% states + constraint are SOS cones
opts.Kx = struct('lin',length(sos1.x));
opts.Kc = struct('sos', length(sos1.g));

% build first solver
S1 = casos.nlsossol('S1','filter-linesearch',sos1,opts);



x0 = [r0;
      x'*x;
      x'*x;
      x'*x;
      x'*x;
      x'*x;
      x'*x];

% call first solver
 sol1 = S1('x0',x0);  

%% plotting
import casos.toolboxes.sosopt.*

Vval0 = sol1.x(1)-sol1.x(end);

figure(2)
% reachability storage function at t = 0s
pcontour(subs(Vval0,x(3:4),zeros(2,1)),0,[-1 1 -2 2],'b')
hold on 
pcontour(subs(rT,x(3:4),zeros(2,1)),0,[-1 1 -2 2],'g--')
% pcontour(subs(r0,x(3:4),zeros(2,1)),0,[-1 1 -2 2]*2,'k')



 % all constraints and decision variables on SDP level
 n_con = S1.stats.single_iterations{end}.conic.size_A.size(1)
 n_dec = S1.stats.single_iterations{end}.conic.size_A.size(1)

