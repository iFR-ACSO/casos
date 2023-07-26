function varargout = degrees(obj)
% Return degrees of polynomial.
%
% Possible syntax:
%
%   deg = degrees(p)
%
% Returns a vector of degrees.
%
%   [mn,mx] = degrees(p)
%
% Returns the minimum and maximum degrees.

if nargout < 2
    varargout = {full(unique(sum(obj.degmat,2)))'};

else
    varargout = {obj.mindeg obj.maxdeg};

end
