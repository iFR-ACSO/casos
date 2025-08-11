function out = str_monoms(obj,omit_constant)
% Return string representation for monomials.

if nargin < 2
    % print constant terms as 1
    omit_constant = false;
end

out = cell(1,obj.nterm);

for it = 1:obj.nterm
    % degrees
    dg = obj.degmat(it,:);

    % monomial
    ms = arrayfun(@exp2str, obj.indets(find(dg)), nonzeros(dg)');
    
    if ~isempty(ms)
        % concatenate exponentials
        out(it) = join(ms,'*');
    elseif ~omit_constant
        % constant
        out{it} = '1';
    end
end

end

function s = exp2str(v,d)
% Convert exponential v^d to string.

    v = str(v);
    switch d
        case 0, s = {''};
        case 1, s = v;
        otherwise, s = compose('%s^%d',v{:},d);
    end
end
