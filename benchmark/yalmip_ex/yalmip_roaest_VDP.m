clear
clc
sdpvar x1 x2  g b

x = [x1;x2];

% system dynamics
f = [-x(2); x(1) + (x(1)^2 - 1)*x(2)];

% Lyapunov function candidate
Vval = 1.5*x(1)^2 - x(1)*x(2) + x(2)^2;
p = x'*x;


l = 1e-6*(x'*x);

% polynomial(indet, maxdeg, mindeg)
[V,cv1]  = polynomial(x,2,2);
[s1,c1]  = polynomial(x,2,2);
[s2,c2]  = polynomial(x,0);

solverset = sdpsettings('solver','mosek', ...
                        'verbose',0, ...
                         'sos.traceobj',0,...   % if not set becomes infeasible quickly; % Minimize trace of Gram matrix in problems without objective function
                         'sos.newton',1,...     % Use Newton polytope to reduce size
                         'sos.congruence',2,... % Block-diagonalize using congruence classes
                         'sos.scale',1);        %scale polynomials




%--------------------------------------------------------------------------
% see https://yalmip.github.io/tutorial/sumofsquaresprogramming/ on how to
% setup constraint sos problems
%--------------------------------------------------------------------------

buse = [];
guse = [];
endTimeParse1 = [];
endTimeParse2 = [];

endTimeSolve1 = [];
endTimeSolve2 = [];
tic
for iter = 1:10

    % find largest stable level set
    lb = -1e2; ub = 1e2;
    % bisection
    startTimeBisec1 = tic;
    counter1 = 0;
    while (ub - lb > 1e-1)
        gtry = (lb+ub)/2;
        
        startTimeParse1 = tic;
        con1 = [sos(s1*(Vval-gtry) - jacobian(Vval,x)*f - l)
                sos(s1)];
        endTimeParse1 = [endTimeParse1 toc(startTimeParse1)];
        
        startTimeSolve1 = tic;
        sol1 = solvesos(con1,[],solverset,c1 );
        endTimeSolve1 = [endTimeSolve1 toc(startTimeSolve1)];

        if sol1.problem == 0
            lb = gtry;
            guse = gtry;
            s1val = replace(s1,c1,value(c1));
        else
            ub = gtry;
        end
           counter1 = counter1 + 1;
    end
    endTimeBisec1(iter) = toc(startTimeBisec1);

         if ~isempty(guse)
            % fprintf('gamma is %g.\n', guse)
        else
             disp(['Problem infeasible in gamma-step in iteration:' num2str(iter)])
            break
        end
   

    % find largest possible shape function
      lb = -1e2; ub = 1e2;
        startTimeBisec2 = tic;
        counter2 = 0;
        while (ub - lb > 1e-1)
        btry = (lb+ub)/2;

        % parse problem
        startTimeParse2 = tic;
        con2 = [sos(s2*(p-btry) +  gtry - Vval)
                sos(s2)];
        endTimeParse2 = [endTimeParse2 toc(startTimeParse2)];
        % solve problem
        startTimeSolve2 = tic;
        sol2 = solvesos(con2,[],solverset,c2 );
        endTimeSolve2 = [endTimeSolve2 toc(startTimeSolve2)];

        if sol2.problem == 0
            lb  = btry;

            % make sure to keep the feasible solution
            buse = btry;
            s2val = replace(s2,c2,value(c2));
            
        else
            ub = btry;
        end
            counter2 = counter2 +1;
        end
        endTimeBisec2(iter) = toc(startTimeBisec2);

        if ~isempty(buse)
            % fprintf('beta is %g.\n', buse)
        else
             disp(['Problem infeasible in beta-step in iteration:' num2str(iter)])
            break
        end

% parse problem
startTimeParse3 = tic;
con3 = [sos(V-l)
        sos(s1val*(V-guse) - jacobian(V,x)*f - l)
        sos(s2val*(p-buse) +  guse - V)
        ];
endTimeParse3(iter) = toc(startTimeParse3);

% solve for Lyapunov function
startTimeSolve3 = tic;
sol3 = solvesos(con3,[],solverset,cv1);
endTimeSolve3(iter) = toc(startTimeSolve3);

if sol3.problem == 0
            Vval = replace(V,cv1,double(cv1));

           fprintf('Iteration %d: b = %g, g = %g.\n',iter,full(buse),full(guse));

            % to make sure we do not use the old solution again
            guse = [];
            buse = [];
else
    disp(['Problem infeasible in V-step in iteration:' num2str(iter)])
    break
end


end
toc

totalTime_parsing = sum(endTimeParse1) + sum(endTimeParse2)
totalTime_solver = sum(endTimeSolve1) + sum(endTimeSolve2) 

totalTime_bisec = sum(endTimeBisec1) + sum(endTimeBisec2)