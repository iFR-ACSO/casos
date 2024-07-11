function r = Fnorm2(obj)
% Compute the squared Frobenius norm of an operator.

r = trace(obj.matrix'*obj.matrix);

end
