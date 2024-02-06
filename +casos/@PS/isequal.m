function tf = isequal(a,b)
% Check if polynomial arrays are equal.

if ~isequal(size(a),size(b))
    tf = false;
    return
end

% else
tf = is_zero(a-b);

end
