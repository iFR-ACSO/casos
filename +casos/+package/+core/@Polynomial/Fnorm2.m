function r = Fnorm2(obj)
% Compute the squared Frobenius norm of a polynomial operator.

assert(is_operator(obj), 'Only allowed for operators.')

r = dot(obj.coeffs,obj.coeffs);

end
