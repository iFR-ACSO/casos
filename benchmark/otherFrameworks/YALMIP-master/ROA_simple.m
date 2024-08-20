clear
clc

sdpvar x1 x2 beta

x = [x1;x2];

x1dot = -x2;
x2dot = x1+(x1^2-1)*x2;
f     = [x1dot; x2dot];


Alin = [0 -1; 1 -1];
Q = eye(2);
P=lyap(Alin',Q);
Vval = x'*P*x;
p = Vval*2;
% Vdot = jacobian(Vval,x)*f;

% sosdecvar always have of polydecvar
[s,c] = polynomial(x,4,2);
[s2,c2] = polynomial(x,4,2);
[V,cv] = polynomial(x,4,1);

solverset = sdpsettings('solver','scs','verbose',0,'sos.traceobj',1);

tic



for i = 1:5


% very simplified bisection

tub = 0;
tlb = -50;
iter = 0;
dopt = [];
while (tub-tlb > 1e-3) || iter >= 100

% Set gamma level for feasibility problem
ttry = (tub+tlb)/2;

sosCon = [ sos(s); 
           sos(Vval - 1e-6*(x'*x)); ...
           sos(-1e-6*(x'*x)+s*(Vval+ttry)-jacobian(Vval,x)*f)];

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

