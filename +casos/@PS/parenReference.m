function varargout = parenReference(obj,indexOp)
% Overwriting matlab.mixin.indexing.RedefinesParen.parenReference

% perform parenthesis reference
idx = indexOp(1);

% select referenced elements
I = logical(sparse(size(obj,1),size(obj,2)));
I.(idx) = true;

if length(indexOp) > 1 && indexOp(2).Type == "Dot"
    % handle getters on referenced polynomial
    done = true;

    switch (indexOp(2).Name)
        case 'mindeg'
            res = min(get_degree(obj,I))';
        case 'maxdeg'
            res = max(get_degree(obj,I))';
        case 'nvars'
            res = length(get_indets(obj,I));
        case 'nterm'
            res = size(get_degmat(obj,I),1);
        case 'indeterminates'
            res = get_indets(obj,I);
        case 'monomials'
            res = get_monoms(obj,I);
        otherwise
            % getter not supported
            done = false;
    end

    if ~done
        % continue
    elseif length(indexOp) > 2
        [varargout{1:nargout}] = res.(indexOp(3:end));
        return
    else
        varargout = {res};
        return
    end
end

% new polynomial
p = casos.PS;

if nnz(I) > 0
    % reference coefficients
    if ~any(cellfun(@(Ind) isnumeric(Ind), idx.Indices))
        % logical or colon reference
        ii = find(I);
    elseif length(idx.Indices) > 1
        % numerical indices
        row = get_indices(size(obj,1),idx.Indices{1});
        col = get_indices(size(obj,2),idx.Indices{2});
        [R,C] = ndgrid(row,col);
        ii = sub2ind(size(obj),R,C);
    else
        ii = get_indices(numel(obj),idx.Indices{1});
    end
    coeffs = obj.coeffs(:,ii);
    
    % remove coefficients, degrees, and/or indeterminates 
    % that do not appear in the referenced polynomial
    [p.coeffs,p.degmat,p.indets] = removeZero(coeffs,obj.degmat,obj.indets);
    % resize
    p.matdim = size(I.(idx));
end

if length(indexOp) > 1
    % forward reference
    [varargout{1:nargout}] = p.(indexOp(2:end));

else
    varargout = {p};
end

end

function ind = get_indices(n,ind)
% Return indices.
    if isequal(ind,':')
        ind = 1:n;
    end
end
