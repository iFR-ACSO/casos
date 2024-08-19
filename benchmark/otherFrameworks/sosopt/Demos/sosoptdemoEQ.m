%---------------------------------------------------------------------
% sosoptdemoEQ
%
% Demonstration of SOSOPT function for solving sos optimizations.
% This example uses SOSOPT to solve a problem with LP, SOS, and 
% equality constraints. 
% See SOSOPT help for more details on the function syntax.
%---------------------------------------------------------------------

% Create polynomial variables, x, and sos decision variable, s(x)
pvar x1 x2;
x = [x1;x2];
[s,D] = sosdecvar('d',x);

% Create SOS and equality constraints on the coefficients of s(x)
D12val = 2.5;
pconstr = [ s>=0; D(1,2)==D12val];

% Add LP constraints on the coefficients of s(x)
LPconstr = [D(1,1)<=4; D(2,2)<=4];
pconstr = [pconstr; LPconstr];

% Create objective function
%th = 2*pi*rand;  
th = 2.4;
c = [cos(th); sin(th)];
obj = c'*[D(1,1); D(2,2)];

% Set options for SOSOPT
opts = sosoptions;
%opts.form = 'image';
%opts.form = 'kernel';
%opts.simplify = 'off';
%opts.solver = 'sdpt3';
%opts.solver = 'csdp';
%opts.solver = 'sdplr';  %NONSYM

% Solve the problem with sosopt
[info,dopt,sossol] = sosopt(pconstr,x,obj,opts);

% Get optimal solution
fprintf('\n----------------Results');
if info.feas
    fprintf('\nSOS/LP Optimization is feasible:\n');
    D11sol = double(subs(D(1,1),dopt))
    D12sol = double(subs(D(1,2),dopt))
    D22sol = double(subs(D(2,2),dopt))
else
    fprintf('\nSOS/LP Optimization is not feasible. \n');    
    return
end

% Plot LP and SOS constraints
figure(1);
plot([4 4 0],[0 4 4],'b'); hold on;
x1dat = linspace(0.5,5);
x2dat = D12val^2./x1dat;
plot(x1dat,x2dat,'r');

% Plot objective function direction
quiver(D11sol,D22sol,c(1),c(2),'k');

% Plot optimal decision variables
plot(D11sol,D22sol,'ko'); hold off;
legend('LP Constraints','SOS+EQ Constraints','Objective Dir.',....
    'Optimal Pt', 'Location','SouthWest')
xlabel('D11');
ylabel('D22');
axis([0.5 4.5 0.5 4.5])
