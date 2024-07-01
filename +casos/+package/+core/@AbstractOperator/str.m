function s = str(obj)
% Return string representation.

sp_in = to_char(obj.sparsity_in);
sp_out = to_char(obj.sparsity_out);

s = compose('Operator: [%s] -> [%s]',sp_in,sp_out);

end
