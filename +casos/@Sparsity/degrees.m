function varargout = degrees(obj,varargin)
% Return degrees of sparsity pattern.
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

dg = get_degree(obj,varargin{:});

if nargout < 2
    varargout = {dg};

else
    varargout = {min(dg) max(dg)};

end

end
