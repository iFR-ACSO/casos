%---------------------------------------------------------------------
% sosoptdemo5
%
% Demonstration of SOSOPT function for solving sos optimizations.
% This example uses SOSOPT to compute an upper bound for the
% structured singular value.
% See SOSOPT help for more details on the function syntax.
%
% Reference:
% S. Prajna, A. Papachristodoulou, P. Seiler, and P.A. Parrilo,
% SOSTOOLS: SOSDEMO5, Section 3.5 of SOSTOOLS Version 2.00
% User's Guide, 2004.
%---------------------------------------------------------------------

% Define polynomial variables
pvar x1 x2 x3 x4 x5 x6 x7 x8;
x = [x1; x2; x3; x4; x5; x6; x7; x8];

% Matrix to be tested
i=sqrt(-1);
alpha = 3 + sqrt(3);
beta = sqrt(3) - 1;
a = sqrt(2/alpha);
b = 1/sqrt(alpha);
c = b;
d = -sqrt(beta/alpha);
f = (1 + i)*sqrt(1/(alpha*beta));
U = [a 0; b b; c i*c; d f];
V = [0 a; b -b; c -i*c; -i*f -d];
M = U*V';

% Constructing A(x)'s
A = cell(4,1);
gam = 0.8724;
Z = monomials(x,1);
for i = 1:4
    H = M(i,:)'*M(i,:) - (gam^2)*sparse(i,i,1,4,4,1);
    H = [real(H) -imag(H); imag(H) real(H)];
    A{i} = (Z.')*H*Z;
end;


% Define SOS decision variables Q1(x),...,Q4(x)
pconstr = polyconstr;
Q = cell(4,1);
for i = 1:4
    Q{i} = sosdecvar(['M' int2str(i)],Z);
    pconstr(i,1) = Q{i}>=0;
end;

% Define constant SOS decision variables rij
R = mpvar('r',4);
idx = find(triu(ones(4),1));
pconstr = [pconstr; R(idx)>=0];

% Define SOS Constraint:
% -sum(Qi(x)*Ai(x)) - sum(rij*Ai(x)*Aj(x)) + I(x) >= 0
expr = 0;
for i = 1:4
    expr = expr - A{i}*Q{i};
end;
for i = 1:4
    for j = (i+1):4
        expr = expr - A{i}*A{j}*R(i,j);
    end;
end;

I = - sum(x.^4);
expr = expr + I;

pconstr(end+1) = expr>=0;

% Set options for SOSOPT
opts = sosoptions;
%opts.form = 'image';
%opts.form = 'kernel';
%opts.solver = 'sdpt3';
%opts.solver = 'csdp';
%opts.solver = 'sdplr';  %NONSYM
%opts.solver = 'sdpam';

% Call sosopt to test feasibility
[info,dopt,sossol] = sosopt(pconstr,x,opts);

% Program is feasible, thus 0.8724 is an upper bound for mu.
fprintf('\n----------------Results');
if info.feas
    fprintf('\nSOS Optimization is feasible:\n');
else
    fprintf('\nSOS Optimization is not feasible. \n');
    return
end

% Verify the Gram Marix Decomposition for last SOS constraint
p = subs(expr,dopt);
z = sossol(end).z;
Q = sossol(end).Q;
e = p - z'*Q*z;
fprintf('\nMax magnitude coefficient of p-z''*Q*z is:')
disp(full(max(abs(e.coefficient))))
fprintf('Minimum eigenvalue of Q is:')
disp(min(eig(Q)));

