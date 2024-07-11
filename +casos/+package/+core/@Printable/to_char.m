function c = to_char(obj)
% Return string representation as character vector.

s = str(obj);

if isscalar(s)
    % return single string representation
    c = s{1};

else
    % no single string representation
    c = sprintf('%dx%d %s',size(s,1),size(s,2),class(obj));
end

end
