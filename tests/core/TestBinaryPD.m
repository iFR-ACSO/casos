% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

classdef (TestTags="PD") TestBinaryPD < TestPolynomialOperations
% Test binary (element-wise) operations on constant polynomials.

properties (SetAccess=protected)
    values       % test polynomials
    references   % reference solutions
end

properties (TestParameter)
    op = {"plus" "minus" "times" "ldivide" "rdivide" "power"};

    dim1 = num2cell(1:6);
    dim2 = num2cell(1:6);
    arg1 = num2cell(1:10);
    arg2 = num2cell(1:10);
end

methods (TestClassSetup)
    function initializeTestData(test_case)
        % Initialize test data for binary operations.
        test_case.loadTestData("binary");

        test_case.fatalAssertLength(test_case.dim1,size(test_case.references.matrix.plus,1));
        test_case.fatalAssertLength(test_case.dim2,size(test_case.references.matrix.plus,2));
        test_case.fatalAssertLength(test_case.arg1,size(test_case.references.scalar.plus,1));
        test_case.fatalAssertLength(test_case.arg2,size(test_case.references.scalar.plus,2));
    end
end

methods (Test, ParameterCombination="pairwise", TestTags="scalar")
    function test_binary_scalar(test_case, op, arg1, arg2)
        % Test binary operation.
        value1 = test_case.values.scalar{1,arg1};
        value2 = test_case.values.scalar{2,arg2};

        reference = test_case.references.scalar.(op){arg1,arg2};

        switch (op)
            case {"plus" "minus" "times"}
                % element-wise addition, subtraction, multiplication
                test_case.evaluate_binary(op,value1,value2,reference);

            case "ldivide"
                % left-divide by constant polynomial
                vars = value1.indeterminates;
                value1_deg0 = 1+subs(value1,vars,ones(length(vars),1));
                
                test_case.evaluate_binary(op,value1_deg0,value2,reference);

            case "rdivide"
                % right-divide by constant polynomial
                vars = value2.indeterminates;
                value2_deg0 = 1+subs(value2,vars,ones(length(vars),1));
                
                test_case.evaluate_binary(op,value1,value2_deg0,reference);

            case "power"
                % element-wise power by degree matrix
                degs = ceil(2*full(project(value2,restrict_terms(sparsity(value2),0))));

                test_case.evaluate_binary(op,value1,degs,reference);

            otherwise
                test_case.assertFail(sprintf("Not implemented: %s.",op));
        end
    end
end

methods (Test, ParameterCombination="pairwise", TestTags=["vector" "inner"])
    function test_binary_inner(test_case, op, dim1, dim2)
        % Test binary operation on row-by-column vectors.
        value1 = test_case.values.vector{1,dim1}';
        value2 = test_case.values.vector{2,dim2};
        
        reference = test_case.references.inner.(op){dim1,dim2};

        switch (op)
            case {"plus" "minus" "times"}
                % element-wise addition, subtraction, multiplication
                test_case.evaluate_binary(op,value1,value2,reference);

            case "ldivide"
                % left-divide by constant polynomial
                vars = value1.indeterminates;
                value1_deg0 = 1+subs(value1,vars,ones(length(vars),1));
                
                test_case.evaluate_binary(op,value1_deg0,value2,reference);

            case "rdivide"
                % right-divide by constant polynomial
                vars = value2.indeterminates;
                value2_deg0 = 1+subs(value2,vars,ones(length(vars),1));
                
                test_case.evaluate_binary(op,value1,value2_deg0,reference);

            case "power"
                % element-wise power by degree matrix
                degs = ceil(2*full(project(value2,restrict_terms(sparsity(value2),0))));

                test_case.evaluate_binary(op,value1,degs,reference);

            otherwise
                test_case.assertFail(sprintf("Not implemented: %s.",op));
        end
    end
end

methods (Test, ParameterCombination="pairwise", TestTags=["vector" "outer"])
    function test_binary_outer(test_case, op, dim1, dim2)
        % Test binary operation on column-by-row vectors.
        value1 = test_case.values.vector{1,dim1};
        value2 = test_case.values.vector{2,dim2}';
        
        reference = test_case.references.outer.(op){dim1,dim2};

        switch (op)
            case {"plus" "minus" "times"}
                % element-wise addition, subtraction, multiplication
                test_case.evaluate_binary(op,value1,value2,reference);

            case "ldivide"
                % left-divide by constant polynomial
                vars = value1.indeterminates;
                value1_deg0 = 1+subs(value1,vars,ones(length(vars),1));
                
                test_case.evaluate_binary(op,value1_deg0,value2,reference);

            case "rdivide"
                % right-divide by constant polynomial
                vars = value2.indeterminates;
                value2_deg0 = 1+subs(value2,vars,ones(length(vars),1));
                
                test_case.evaluate_binary(op,value1,value2_deg0,reference);

            case "power"
                % element-wise power by degree matrix
                degs = ceil(2*full(project(value2,restrict_terms(sparsity(value2),0))));

                test_case.evaluate_binary(op,value1,degs,reference);

            otherwise
                test_case.assertFail(sprintf("Not implemented: %s.",op));
        end
    end
end

methods (Test, ParameterCombination="pairwise", TestTags=["matrix"])
    function test_binary_matrix(test_case, op, dim1, dim2)
        % Test binary operation on matrices.
        value1 = test_case.values.matrix{1,dim1};
        value2 = test_case.values.matrix{2,dim2};

        reference = test_case.references.matrix.(op){dim1,dim2};

        switch (op)
            case {"plus" "minus" "times"}
                % element-wise addition, subtraction, multiplication
                test_case.evaluate_binary(op,value1,value2,reference);

            case "ldivide"
                % left-divide by constant polynomial
                vars = value1.indeterminates;
                value1_deg0 = 1+subs(value1,vars,ones(length(vars),1));
                
                test_case.evaluate_binary(op,value1_deg0,value2,reference);

            case "rdivide"
                % right-divide by constant polynomial
                vars = value2.indeterminates;
                value2_deg0 = 1+subs(value2,vars,ones(length(vars),1));
                
                test_case.evaluate_binary(op,value1,value2_deg0,reference);

            case "power"
                % element-wise power by degree matrix
                degs = ceil(2*full(project(value2,restrict_terms(sparsity(value2),0))));

                test_case.evaluate_binary(op,value1,degs,reference);

            otherwise
                test_case.assertFail(sprintf("Not implemented: %s.",op));
        end
    end
end

methods (Test, ParameterCombination="pairwise", TestTags=["monomials"])
    function test_binary_monomials(test_case, op, dim1, dim2)
        % Test binary operation on monomial matrices.
        value1 = test_case.values.matrix{4,dim1};
        value2 = test_case.values.matrix{2,dim2};

        reference = test_case.references.monomials.(op){dim1,dim2};

        switch (op)
            case {"plus" "minus" "times"}
                % element-wise addition, subtraction, multiplication
                test_case.evaluate_binary(op,value1,value2,reference);

            case "ldivide"
                % left-divide by constant polynomial
                vars = value1.indeterminates;
                value1_deg0 = 1+subs(value1,vars,ones(length(vars),1));
                
                test_case.evaluate_binary(op,value1_deg0,value2,reference);

            case "rdivide"
                % right-divide by constant polynomial
                vars = value2.indeterminates;
                value2_deg0 = 1+subs(value2,vars,ones(length(vars),1));
                
                test_case.evaluate_binary(op,value1,value2_deg0,reference);

            case "power"
                % element-wise power by degree matrix
                degs = ceil(2*full(project(value2,restrict_terms(sparsity(value2),0))));

                test_case.evaluate_binary(op,value1,degs,reference);

            otherwise
                test_case.assertFail(sprintf("Not implemented: %s.",op));
        end
    end
end

methods
    function evaluate_binary(test_case, op, value1, value2, reference)
        % Evaluate binary operation.
        if (~isrow(value1) && ~isrow(value2) && size(value1,1) ~= size(value2,1)) ...
                    || (~iscolumn(value1) && ~iscolumn(value2) && size(value1,2) ~= size(value2,2))
            % dimension mismatch
            diagtext = sprintf('Dimension mismatch expected: %s vs. %s.',mat2str(size(value1)),mat2str(size(value2)));
            test_case.verifyError(@() feval(op,value1,value2),?casos.package.core.IncompatibleSizesError,diagtext);
            return
        end

        % else
        actual = feval(op,value1,value2);

        % perform assertion
        test_case.verifyEqualPolynomial(actual,reference,"RelTol",1e-14);
    end
end

end
