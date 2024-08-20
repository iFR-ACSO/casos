%---------------------------------------------------------------------
% sosoptdemoEQ
%
% Demonstration of SOSOPT function for solving sos optimizations. This
% example uses SOSOPT to solve a problem with LP and SOS constraints.
% See SOSOPT help for more details on the function syntax.
%---------------------------------------------------------------------

% Create polynomial variables, x, and decision variables, d
pvar d1 d2 x1 x2;
x = [x1;x2];
d = [d1;d2];

% Create an SOS constraint on the coefficients of s(x)
s = d1*x1^2+4*x1*x2+d2*x2^2;
pconstr = s>=0;

% Create an LP constraint on the coefficients
if 1
    % Group into one LP constraint 
    A = [-1 0; 1 0; 0 -1; 0 1];
    b = [-1; 4; -1; 4];
    LPconstr = A*d<=b;
else
    % Specify LP constraints one at a time and stack vertically
    LPconstr = [ d1>=1; d1<=4; d2>=1; d2<=4];
end
pconstr = [pconstr; LPconstr];

% Create objective function
th = 2*pi*rand;  
%th = 5*pi/4; 
c = [cos(th); sin(th)];
obj = c'*d;

% Set options for SOSOPT
opts = sosoptions;
%opts.form = 'image';
opts.form = 'kernel';
%opts.solver = 'sdpt3';
%opts.solver = 'csdp';
%opts.solver = 'sdplr';

% Solve the problem with sosopt
[info,dopt,sossol] = sosopt(pconstr,x,obj,opts);

% Get optimal solution
fprintf('\n----------------Results');
if info.feas
    fprintf('\nSOS/LP Optimization is feasible:\n');
    d1sol = double(subs(d1,dopt))
    d2sol = double(subs(d2,dopt))
else
    fprintf('\nSOS/LP Optimization is not feasible. \n');    
    return
end

% Plot LP and SOS constraints
figure(1);
plot([1 4 4 1 1],[1 1 4 4 1],'b'); hold on;
x1dat = linspace(0.5,5);
x2dat = 4./x1dat;
plot(x1dat,x2dat,'r');

% Plot objective function direction
quiver(d1sol,d2sol,c(1),c(2),'k');

% Plot optimal decision variables
plot(d1sol,d2sol,'ko'); hold off;
legend('LP Constraints','SOS Constraint','Objective Dir.',....
    'Optimal Pt', 'Location','SouthWest')
xlabel('d1');
ylabel('d2');
axis([0.5 4.5 0.5 4.5])
