function argout = eval(obj,argin)
% Call SeDuMi interface.

% pre-process bounds
lba = sparse(argin{4});
uba = sparse(argin{5});
cba = sparse(argin{6});
lbx = sparse(argin{7});
ubx = sparse(argin{8});
cbx = sparse(argin{9});

% dimensions of original problem
nl = length(lbx);
nc = length(cbx);
ml = length(lba);
mc = length(cba);

% detect infinite lower variable bounds
If = find(isinf(lbx));
% prepare for infinite lower variable bounds
argin{7}(If) = 0;

% evaluate problem structure
prob = call(obj.fhan,argin);

% to double
A = sparse(prob{1});
b = sparse(prob{2});
c = sparse(prob{3});
% cone
K = obj.cone;

% options to SeDuMi
opts = obj.sdpopt.solveroptions;

% reorder decision variables
idx = [If' setdiff(1:length(c),If)];

% remove trivial constraints
I = false(size(b));
J = false(size(c));

% detect equality constraints
Ila = find(lba == uba);
% remove lower bound constraints
I(Ila) = true;
% remove slack variables (sua,sla)
J(nl+[Ila; ml+Ila]) = true;

% detect constant variables
Ilx = find(lbx == ubx);
% remove slack variable sux
% (this leaves zl(i) = 0)
J(nl+2*ml+Ilx) = true;

% detect infinite bounds
Iba = find(isinf([uba;lba]));
Ibx = find(isinf(ubx));
% remove infinite bound constraints
I(Iba) = true;
I(2*ml+mc+Ibx) = true;
% remove slack variables (sua,sla,sux)
J(nl+[Iba; 2*ml+Ibx]) = true;

% purge constraints
A(I,:) = [];
b(I)   = [];

% purge variables
idx(J) = [];

A = A(:,idx);
c = c(idx);

% modify cone
K.f = length(If);
K.l = K.l - nnz(J) - length(If);

% call SeDuMi
[x_,y_,info] = sedumi(A,b,full(c),K,opts);

% assign full solution
x = sparse(idx,1,x_,length(J),1);
y = sparse(find(~I),1,y_,length(I),1);

if ~obj.sdpopt.error_on_fail
    % continue regardless of feasibility
elseif info.pinf
    % primal infeasible
    error('Conic problem is primal infeasible.')
elseif info.dinf
    % dual infeasible
    error('Conic problem is dual infeasible.')
elseif info.numerr
    % numerical errors
    error('Optimizer run into numerical error (numerr=%d, feasratio=%g)',info.numerr,info.feasratio)
end

% parse solution
argout = call(obj.ghan,[argin {x y}]);

end