% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

classdef (TestTags="PD") TestPoly2BasisPD < TestPolynomialOperations
% Test poly2basis operation on constant polynomials.

properties (SetAccess=protected)
    values       % test polynomials
    references   % reference solutions
end

properties (TestParameter)
    dim1 = num2cell(1:6);
    dim2 = num2cell(1:6);
    arg1 = num2cell(1:10);
    arg2 = num2cell(1:10);
end

methods (TestClassSetup)
    function initializeTestData(test_case)
        % Initialize test data for poly2basis operation.
        test_case.loadTestData("poly2basis");

        test_case.fatalAssertLength(test_case.dim1,size(test_case.references.matrix,1));
        test_case.fatalAssertLength(test_case.dim2,size(test_case.references.matrix,2));
        test_case.fatalAssertLength(test_case.arg1,size(test_case.references.scalar,1));
        test_case.fatalAssertLength(test_case.arg2,size(test_case.references.scalar,2));
    end
end

methods (Test, ParameterCombination="pairwise", TestTags="scalar")
    function test_poly2basis_scalar(test_case, arg1, arg2)
        % Test poly2basis operation.
        value = test_case.values.scalar{1,arg1};
        basis = sparsity(densify(test_case.values.scalar{2,arg2})); % project onto dense basis

        reference = test_case.references.scalar{arg1,arg2};

        test_case.evaluate_poly2basis(value,basis,reference);
    end
end

methods (Test, ParameterCombination="pairwise", TestTags=["vector" "column"])
    function test_poly2basis_column(test_case, dim1, dim2)
        % Test poly2basis operation on column vectors.
        value = test_case.values.vector{1,dim1};
        basis = sparsity(densify(test_case.values.vector{2,dim2})); % project onto dense basis
        
        reference = test_case.references.column{dim1,dim2};

        test_case.evaluate_poly2basis(value,basis,reference);
    end
end

methods (Test, ParameterCombination="pairwise", TestTags=["vector" "row"])
    function test_poly2basis_row(test_case, dim1, dim2)
        % Test poly2basis operation on row vectors.
        value = test_case.values.vector{1,dim1}';
        basis = sparsity(densify(test_case.values.vector{2,dim2}')); % project onto dense basis
        
        reference = test_case.references.row{dim1,dim2};

        test_case.evaluate_poly2basis(value,basis,reference);
    end
end

methods (Test, ParameterCombination="pairwise", TestTags="matrix")
    function test_poly2basis_matrix(test_case, dim1, dim2)
        % Test poly2basis operation on matrices.
        value = test_case.values.matrix{1,dim1};
        basis = sparsity(densify(test_case.values.matrix{2,dim2})); % project onto dense basis
        
        reference = test_case.references.matrix{dim1,dim2};

        test_case.evaluate_poly2basis(value,basis,reference);
    end
end

methods
    function evaluate_poly2basis(test_case, value, basis, reference)
        % Evaluate poly2basis operation.
        if ~isscalar(value) && ~(isequal(size(value),[0 0]) && isempty(basis)) && ~isequal(size(value),size(basis))
            % size mismatch
            diagtext = sprintf('Dimension mismatch expected: %s vs. %s.',mat2str(size(value)),mat2str(size(basis)));
            test_case.verifyError(@() poly2basis(value,basis),?casos.package.core.IncompatibleSizesError,diagtext);
            return
        end

        % else
        actual = poly2basis(value,basis);

        % perform assertion
        test_case.verifyClass(actual,?casadi.DM);
        test_case.verifyEqual(full(actual),reference,"RelTol",1e-15);
    end
end

end
