function Lz = newton_reduce(obj, Pdegmat, Zdegmat)
% Removes monomials outside half Newton polytope
% Strategy inpired from: Simplification Methods for Sum-of-Squares 
% Programs, Peter Seiler et al.

% options for linprog
options = optimoptions('linprog');
options.Display = 'off';

% build LP (part 1)
bfixed = [zeros(size(Pdegmat,1),1); 1];

% trivial removal of monomials
keep_trivial = ismember(Zdegmat*2,Pdegmat,'rows');
keep = true(size(keep_trivial));

% try to go over each possible monomial basis and verify if it belongs to
% the newton polytope by checking for the existance of a hyperplane 
for i = 1:length(keep_trivial)

    if keep_trivial(i) || ~keep(i)
        continue;
    end
    
    % build LP (part 2)
    q = Zdegmat(i,:)*2;
    c = ([-q 1])';
    F_struc = ([bfixed -[Pdegmat -ones(size(Pdegmat,1),1); q -1]]);
    
    A = -F_struc(1:end,2:end);
    b = F_struc(1:end,1);   

    % solve LP
    [x,~,flag] = linprog(c, A, b, [], [], [], [], options);
    
    % in case the LP is not feasible, giving an empty output 'x' 
    if isempty(x)
        x = zeros(length(c),1);
    end
    
    % If the LP gives an unbounded solution or the only solution is a
    % zero vector
    if (flag > 0 && ([-q 1]*x(:) < 0)) || flag == -3
        a = x(1:end-1);
        b = x(end);
        u = 2*Zdegmat*a - b > sqrt(eps);
        keep(u) = 0;
    end
end

Lz = keep;
Lz = sparse(Lz);
end
