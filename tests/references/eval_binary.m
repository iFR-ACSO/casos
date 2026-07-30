% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function res = eval_binary(op,arg1,arg2)
% Evaluate binary operation.

if isrow(arg1)
    % repeat to match rows of second argument
    arg1 = eval_repmat(arg1,size(arg2,1),1);
elseif isrow(arg2)
    % repeat to match rows of first argument
    arg2 = eval_repmat(arg2,size(arg1,1),1);
else
    assert(size(arg1,1) == size(arg2,1), 'Dimension mismatch.')
end

if iscolumn(arg1)
    % repeat to match columns of second argument
    arg1 = eval_repmat(arg1,1,size(arg2,2));
elseif iscolumn(arg2)
    % repeat to match columns of first argument
    arg2 = eval_repmat(arg2,1,size(arg1,2));
else
    assert(size(arg1,2) == size(arg2,2), 'Dimension mismatch.')
end

switch (op)
    case {"plus" "minus" "times"}
        % element-wise addition, subtraction, multiplication
        res = feval(op,arg1,arg2);

    case "ldivide"
        % left-side division B.\A
        vars = arg1.varname;
        arg1_double = 1+double(subs(arg1,vars,ones(size(vars))));
        res = times((arg1_double).\1,arg2);

    case "rdivide"
        % right-side division A./B
        vars = arg2.varname;
        arg2_double = 1+double(subs(arg2,vars,ones(size(vars))));
        res = rdivide(arg1,arg2_double);

    case "power"
        % binary power
        degs = ceil(2*double(cleanpoly(arg2,[],0)));
        res = power(arg1,degs);

    otherwise
        error('Not implemented: %s.', op)
end

end
