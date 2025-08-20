function argout = eval_on_basis(obj,argin)
% Base implementation of simple sequential

%import casos.package.UnifiedReturnStatus

args = argin;

xi_k = args{1};

% parameter from nonlinear problem
p0   = args{2};


Bk = eye(obj.sizeHessian(1));

% initialize iteration info
info = cell(1,obj.opts.max_iter);


args{2}  = [p0; xi_k; Bk(:);zeros(length(xi_k),1)];
% evaluate convex SOS problem
sol = eval_on_basis(obj.sossolver, args);

% store iteration info
info{i} = obj.sossolver.get_stats;


% return last solution
argout = sol;

% store iteration info
obj.info.iter = info;

end
