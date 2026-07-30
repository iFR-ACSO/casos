% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

classdef (TestTags="PD") TestRemoveCoeffs < TestPolynomialOperations
% Test remove_coeffs operation on constant polynomials.

properties (SetAccess=protected)
    values       % test polynomials
    references   % reference solutions
end

properties (TestParameter)
    decade = num2cell(1:3);
    dim = num2cell(1:6);
    arg = num2cell(1:10);
end

methods (TestClassSetup)
    function initializeTestData(test_case)
        % Initialize test data for remove_coeff operation.
        test_case.loadTestData("remove_coeffs");

        test_case.fatalAssertLength(test_case.decade,size(test_case.references.scalar.remove_coeffs,2));
        test_case.fatalAssertLength(test_case.dim,size(test_case.references.matrix.remove_coeffs,1));
        test_case.fatalAssertLength(test_case.arg,size(test_case.references.scalar.remove_coeffs,1));
    end
end

methods (Test, ParameterCombination="pairwise", TestTags="scalar")
    function test_remove_coeffs_scalar(test_case, decade, arg)
        % Test remove_coeffs operation.
        value1 = test_case.values.scalar{1,arg};
        value2 = test_case.values.scalar{2,arg};

        reference = test_case.references.scalar.remove_coeffs{arg,decade};

        test_case.evaluate_remove_coeffs(decade,value1,value2,reference);
    end
end

methods (Test, ParameterCombination="pairwise", TestTags=["vector" "column"])
    function test_remove_coeffs_column(test_case, decade, dim)
        % Test remove_coeffs operation on column vectors.
        value1 = test_case.values.vector{1,dim};
        value2 = test_case.values.vector{2,dim};

        reference = test_case.references.column.remove_coeffs{dim,decade};

        test_case.evaluate_remove_coeffs(decade,value1,value2,reference);
    end
end

methods (Test, ParameterCombination="pairwise", TestTags=["vector" "row"])
    function test_remove_coeffs_row(test_case, decade, dim)
        % Test remove_coeffs operation on row vectors.
        value1 = test_case.values.vector{1,dim}';
        value2 = test_case.values.vector{2,dim}';

        reference = test_case.references.row.remove_coeffs{dim,decade};

        test_case.evaluate_remove_coeffs(decade,value1,value2,reference);
    end
end

methods (Test, ParameterCombination="pairwise", TestTags="matrix")
    function test_remove_coeffs_matrix(test_case, decade, dim)
        % Test remove_coeffs operation on matrices.
        value1 = test_case.values.matrix{1,dim};
        value2 = test_case.values.matrix{2,dim};

        reference = test_case.references.matrix.remove_coeffs{dim,decade};

        test_case.evaluate_remove_coeffs(decade,value1,value2,reference);
    end
end

methods
    function evaluate_remove_coeffs(test_case, decade, value1, value2, reference)
        % Evaluate remove_coeffs operation.
        tol = prctile([0; full(poly2basis(densify(value2)))],10*decade); % project onto dense basis
        actual = remove_coeffs(value1,tol);

        % perform assertion
        test_case.verifyEqualPolynomial(actual,reference,"RelTol",1e-12);
    end
end

end
