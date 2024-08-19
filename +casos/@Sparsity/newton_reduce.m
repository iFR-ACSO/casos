function Lz = newton_reduce(obj, Pdegmat, Zdegmat)
% removes monomials outside half Newton polytope

% options for linprog
options = optimoptions('linprog');
options.Display = 'off';

% build LP (part 1)
bfixed = [zeros(size(Pdegmat,1),1); 1];

% trivial removal of monomials
keep = ismember(Zdegmat*2,Pdegmat,'rows');
dontkeep = false(size(keep));

% try to go over each possible monomial basis and verify of it belongs to
% the newton polytope by checking for the existance of a hyperplane 
for i = 1:length(keep)
    
    if ~keep(i) && ~dontkeep(i)
        % build LP (part 2)
        q = Zdegmat(i,:)'*2;
        c = ([-q' 1])';
        F_struc = ([bfixed -[Pdegmat -ones(size(Pdegmat,1),1); q' -1]]);
        
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
        if (flag > 0 && ([-q' 1]*x(:) < 0)) || flag == -3
            a = x(1:end-1);
            b = x(end);
            u = find(a'*2*Zdegmat' - b > sqrt(eps));
            dontkeep(u) = 1;
        end
    end  
end

Lz = ~dontkeep;
Lz = sparse(Lz);
end
