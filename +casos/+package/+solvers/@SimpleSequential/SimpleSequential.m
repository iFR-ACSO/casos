classdef SimpleSequential < casos.package.solvers.SequentialCommon
% A simple sequential sum-of-squares algorithm.

properties (SetAccess=private)
    class_name = 'SimpleSequential';
end


methods (Access=protected)
    % iteration for overloading
    varargout = run_iteration(varargin);
    varargout = do_single_iteration(varargin);

    % internal evaluation
    argout = eval_on_basis(obj,argin);

end


end
