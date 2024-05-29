function s = str(obj)
% Return string representation.

sp_in = signature(obj.sparsity_in);
sp_out = signature(obj.sparsity_out);

s = compose('Operator: [%s] -> [%s]',sp_in{:},sp_out{:});

end
