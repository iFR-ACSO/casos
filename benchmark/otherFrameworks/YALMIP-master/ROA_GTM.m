clear
clc

sdpvar x1 x2 x3 x4  u1 u2

x = [x1;x2;x3;x4];
u = [u1; u2];

gtm_dyn

% sdisplay(f(1))
% sdisplay(f(2))
% sdisplay(f(3))
% sdisplay(f(4))
% Vdot = jacobian(Vval,x)*f;

% sosdecvar always have of polydecvar
[s,c] = polynomial(x,1);
[s2,c2] = polynomial(x,4,2);
[V,cv] = polynomial(x,4,1);

solverset = sdpsettings('solver','mosek','verbose',0,'sos.clean',1e-6);

tic



for i = 1:5


% very simplified bisection

tub = 0;
tlb = -1000;
iter = 0;
dopt = [];
while (tub-tlb > 1e-3) || iter >= 1000

% Set gamma level for feasibility problem
ttry = (tub+tlb)/2;
% ttry = -0.186920;
sosCon = [ sos(s); 
           sos(s*(Vval+ttry)-jacobian(Vval,x)*f -1e-6*(x'*x))];

[sol,m,Q,residuals,everything] = solvesos(sosCon,...
                                            [],...
                                            solverset,...
                                            c);
    % Update bisection bounds
    if sol.problem == 0
        tub = ttry;
        dopt = ttry;
        sval = replace(s,c,double(c));
        solvertime = sol.solvertime;

    else
        tlb = ttry;
    end
    
    iter = iter+1;

end

if ~isempty(dopt)
gamma = -dopt

else
     disp(['gamma step infeasible in iteration ' num2str(i)])
    break
end



tub = 0;
tlb = -50;
iter = 0;
dopt = [];
% very simplified bisection
while (tub-tlb > 1e-3) || iter >= 100

% Set gamma level for feasibility problem
ttry = (tub+tlb)/2;

sosCon = [ sos(s2); 
           sos(s2*(p+ttry)-(Vval-gamma))];

[sol,m,Q,residuals,everything] = solvesos(sosCon,...
                                            [],...
                                            solverset,...
                                            c2);
    % Update bisection bounds
    if sol.problem == 0
        tub = ttry;
        dopt = ttry;
        s2val = replace(s2,c2,double(c2));
        solvertime = sol.solvertime;

    else
        tlb = ttry;
    end
    
    iter = iter+1;

end


if ~isempty(dopt)
beta = -dopt

else
     disp(['Beta step infeasible in iteration ' num2str(i)])
    break
end


sosCon = [ sos(V - 1e-6*(x'*x)); ...
          sos(s2val*(p-beta)-(V-gamma))
           sos(-1e-6*(x'*x)+sval*(V-gamma)-jacobian(V,x)*f)];

[sol,m,Q,residuals,everything] = solvesos(sosCon,...
                                            [],...
                                            solverset,...
                                            cv);

if sol.problem == 0
Vval = replace(V,cv,double(cv));
else
    disp(['V step infeasible in iteration ' num2str(i)])
    break
end
end


totaltime = toc

