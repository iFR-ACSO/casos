function [x, p, fval] = findtrim(f, X0, P0, cond, opt, bnds)
% Finds trim condition under constraints.
%
%% Usage and description
%
%   [x*, p*] = findtrim(f, x0, p0, cond)
%
% Finds trim condition (x*,p*) of the non-linear function |f(x,p)| under
% the constraints |cond(x,p) == 0| starting at the initial guess (x0,p0).
%
%   [x*, p*] = findtrim(f, x0, p0, (cond | []), opt, [bnds])
%
% Additionally solves the optimization problem (x*,p*) = min opt(x,p) s.t.
% (x,p) is subject to trim condition, optional constraints, and/or given
% bounds |bnds|. Bounds must be a k-by-2 matrix where the left column
% represents lower bound, right column upper bounds. 
% Set bnds_(i,j) to +/- Inf if unbounded.
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      mailto:torbjoern.cunis@onera.fr
% * Created:    2018-05-17
% * Changed:    2018-05-17
%
%%

if ~exist('cond', 'var')
    cond = @(X,P) P - P0;
elseif isempty(cond)
    cond = @(~,~) [];
end

if ~exist('opt', 'var')
    opt = [];
elseif ~exist('bnds', 'var')
    bnds = zeros(0,2);
else
    assert(size(bnds,2) == 2, 'Bounds must be k-by-2 matrix [lb ub].');
end

% number of states
n = length(X0);
% number of parameters
m = length(P0);

% states
X = @(xp) xp(1:n);
% parameters
P = @(xp) xp(n+(1:m));

if isempty(opt)
    % zero problem structure
    problem.objective = @(xp) [cond(X(xp), P(xp)); f(X(xp), P(xp))];
    problem.x0 = [X0; P0];
    problem.solver  = 'fsolve';
    problem.options = optimoptions('fsolve', 'Algorithm', 'trust-region-reflective',...
                                   'FunctionTolerance',1e-20,'StepTolerance',1e-20,...
                                   'MaxFunctionEvaluations',3000,...
                                   'OptimalityTolerance',1e-17);

    % find trim condition for f(x,p) == 0 & g(x,p) == 0
    [xp0,fval] = fsolve(problem);
else
    % optimization problem structure
    problem.objective = @(xp) opt(X(xp), P(xp));
    problem.x0 = [X0; P0];
    problem.lb = bnds(:,1);
    problem.ub = bnds(:,2);
    problem.nonlcon = @eqcon;
    problem.solver  = 'fmincon';
    problem.options = optimoptions('fmincon');
    
    % find optimal trim condition
    [xp0,fval] = fmincon(problem);
end

x = X(xp0);
p = P(xp0);


function [ciq,ceq] = eqcon(xp)
% Nonlinear constraints for fmincon.

    ciq = [];
    ceq = [cond(X(xp), P(xp)); f(X(xp), P(xp))];
end

end