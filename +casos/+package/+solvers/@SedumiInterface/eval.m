function argout = eval(obj,argin)
% Call SeDuMi interface.

prob = call(obj.fhan,argin);

% to double
A = sparse(prob{1});
b = sparse(prob{2});
c = sparse(prob{3});
% cone
K = obj.cone;

% options to SeDuMi
opts = obj.sdpopt.solveroptions;

% detect infinite bounds
I = isinf(b);
% remove infinite interval constraints
A(I,:) = [];
b(I)   = [];

% pre-init dual solution
y = spalloc(length(I),1,nnz(~I));

% call SeDuMi
[x,y(~I),info] = sedumi(A,b,c,K,opts);

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