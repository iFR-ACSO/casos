% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

classdef (TestTags="PS") TestDotPS < TestSymbolicOperations
% Test dot product on symbolic polynomials.

properties (SetAccess=protected)
    values       % test polynomials
    references   % reference solutions
end

properties (TestParameter)
    symb1 = {false true};
    symb2 = {false true};

    dim1 = num2cell(1:6);
    dim2 = num2cell(1:6);
    arg1 = num2cell(1:10);
    arg2 = num2cell(1:10);
end

methods (TestClassSetup)
    function initializeTestData(test_case)
        % Initialize test data for dot product.
        test_case.loadTestData("dot");

        test_case.fatalAssertLength(test_case.dim1,size(test_case.references.matrix.dot,1));
        test_case.fatalAssertLength(test_case.dim2,size(test_case.references.matrix.dot,2));
        test_case.fatalAssertLength(test_case.arg1,size(test_case.references.scalar.dot,1));
        test_case.fatalAssertLength(test_case.arg2,size(test_case.references.scalar.dot,2));
    end
end

methods (Test, ParameterCombination="pairwise", TestTags="scalar")
    function test_dot_scalar(test_case, symb1, symb2, arg1, arg2)
        % Test dot product on scalars.
        value1 = test_case.values.scalar{1,arg1};
        value2 = test_case.values.scalar{2,arg2};

        reference = test_case.references.scalar.dot{arg1,arg2};

        test_case.evaluate_dot(symb1,symb2,value1,value2,reference);
    end
end

methods (Test, ParameterCombination="pairwise", TestTags="matrix")
    function test_dot_matrix(test_case, symb1, symb2, dim1, dim2)
        % Test dot product on matrices.
        value1 = test_case.values.matrix{1,dim1};
        value2 = test_case.values.matrix{2,dim2};

        reference = test_case.references.matrix.dot{dim1,dim2};

        test_case.evaluate_dot(symb1,symb2,value1,value2,reference);
    end
end

methods
    function evaluate_dot(test_case, symb1, symb2, value1, value2, reference)
        % Evaluate dot product.
        test_case.assumeTrue(symb1 || symb2, "Operands not symbolic.");
        
        % symbolic polynomials
        [p1,symbol1,argument1] = test_case.get_operand(symb1,value1);
        [p2,symbol2,argument2] = test_case.get_operand(symb2,value2);

        if ~isequal(size(value1), size(value2))
            % size mismatch
            diagtext = sprintf('Dimension mismatch expected: %s vs. %s.',mat2str(size(value1)),mat2str(size(value2)));
            test_case.verifyError(@() dot(p1,p2),?casos.package.core.IncompatibleSizesError,diagtext);
            return
        end

        % else:
        % build symbolic function
        expression = dot(p1,p2);
        f = casos.Function('f',[symbol1 symbol2],{expression});

        actual = f(argument1{:},argument2{:});

        % perform assertion
        test_case.verifyClass(expression,?casadi.SX);
        test_case.verifyEqual(full(actual),reference,"RelTol",1e-15);
    end
end

end
