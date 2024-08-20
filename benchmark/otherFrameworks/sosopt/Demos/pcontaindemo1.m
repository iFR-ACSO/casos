%---------------------------------------------------------------------
% pcontaindemo1
%
% Demonstration of PCONTAIN function for proving set containments
% using SOS optimizations.  The example uses PCONTAIN to maximize
% the size of a circle inside the contour of a 6th degree polynomial.
% See PCONTAIN help for more details on the function syntax.
%
% Reference: 
% P. Seiler and G. Balas, Quasiconvex Sum-of-Squares Programming,
% IEEE Conference on Decision and Control, 2010.
%---------------------------------------------------------------------

% Create polynomial variables
pvar x1 x2;
x = [x1;x2];

% Construct poly for S1 := { x : g1(x)<= 0}
g1 = 0.3*x1^6 + 0.05*x2^6 - 0.5*x1^5 - 1.4*x1^3*x2 ...
    + 2.3*x1^2*x2^2 - 0.9*x1^3 + 2.6*x1^2*x2 - 1;

% Construct poly for S2 := { x : g2(x)<= gamma}
g2 = x'*x;

% Define monomials for S-procedure multiplier, s(x)
z = monomials(x,0:2);

% Set pcontain options
opts = gsosoptions;
opts.minobj = 0;
opts.maxobj = 1000;

% Use pcontain to find the largest gamma such that S2 is a subset of S1.
% gbnds gives lower/upper bounds on optimal gamma. sopt is the optimal
% S-procedure multiplier.
[gbnds,sopt] = pcontain(g1,g2,z,opts)
gamma = gbnds(1);

% Plot contours of maximal disk and set S1.
plotdomain = [-2 3 -2 2];
pcontour(g1,0,plotdomain,'b'); 
hold on;
pcontour(g2,gamma,plotdomain,'r')
axis equal; axis(plotdomain)

