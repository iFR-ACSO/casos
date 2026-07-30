% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

classdef (TestTags="PS") TestMtimesPS < TestSymbolicOperations
% Test matrix multiplication operations on symbolic polynomials.

properties (SetAccess=protected)
    values       % test polynomials
    references   % reference solutions
end

properties (TestParameter)
    op = {"mtimes" "kron" "mldivide" "mrdivide"};
    symb1 = {false true};
    symb2 = {false true};
    pow = num2cell(0:3);

    dim1 = num2cell(1:6);
    dim2 = num2cell(1:6);
    arg1 = num2cell(1:10);
    arg2 = num2cell(1:10);
end

methods (TestClassSetup)
    function initializeTestData(test_case)
        % Initialize test data for matrix multiplication operations.
        test_case.loadTestData("mtimes");

        test_case.fatalAssertLength(test_case.pow,size(test_case.references.scalar.mpower,2));
        test_case.fatalAssertLength(test_case.dim1,size(test_case.references.matrix.mtimes,1));
        test_case.fatalAssertLength(test_case.dim2,size(test_case.references.matrix.mtimes,2));
        test_case.fatalAssertLength(test_case.arg1,size(test_case.references.scalar.mtimes,1));
        test_case.fatalAssertLength(test_case.arg2,size(test_case.references.scalar.mtimes,2));
    end
end

methods (Test, ParameterCombination="pairwise", TestTags="scalar")
    function test_mtimes_scalar(test_case, op, symb1, symb2, arg1, arg2)
        % Test matrix operation on scalars.
        value1 = test_case.values.scalar{1,arg1};
        value2 = test_case.values.scalar{2,arg2};
        
        reference = test_case.references.scalar.(op){arg1,arg2};

        switch (op)
            case {"mtimes" "kron"}
                % matrix multiplication, Kronecker product
                test_case.evaluate_mtimes(op,symb1,symb2,value1,value2,reference);

            case {"mldivide"}
                % matrix left-divide
                vars = value1.indeterminates;
                value1_deg0 = 1+subs(value1,vars,ones(length(vars),1));
        
                test_case.evaluate_mtimes(op,symb1,symb2,value1_deg0,value2,reference);

            case {"mrdivide"}
                % matrix right-divide
                vars = value2.indeterminates;
                value2_deg0 = 1+subs(value2,vars,ones(length(vars),1));
        
                test_case.evaluate_mtimes(op,symb1,symb2,value1,value2_deg0,reference);
        end
    end

    function test_mpower_scalar(test_case, pow, arg1)
        % Test matrix power operation on scalars (scalar exponent).
        value = test_case.values.scalar{1,arg1};
        
        % symbolic polynomial
        [p,symbol,argument] = test_case.get_operand(true,value);

        % build symbolic function
        expression = mpower(p,pow);
        f = casos.Function('f',symbol,{expression});

        actual = f(argument{:});
        reference = test_case.references.scalar.mpower{arg1,pow+1};

        % perform assertion
        test_case.verifyEqualPolynomial(actual,reference,"RelTol",1e-12)
    end
end

methods (Test, ParameterCombination="pairwise", TestTags=["vector" "inner"])
    function test_mtimes_inner(test_case, symb1, symb2, dim1, dim2)
        % Test matrix multiplication on row-by-column vectors.
        value1 = test_case.values.vector{1,dim1}';
        value2 = test_case.values.vector{2,dim2};

        reference = test_case.references.inner.mtimes{dim1,dim2};

        test_case.evaluate_mtimes("mtimes",symb1,symb2,value1,value2,reference);
    end

    function test_kron_inner(test_case, symb1, symb2, dim1, dim2)
        % Test Kronecker product on row-by-column vectors.
        value1 = test_case.values.vector{1,dim1}';
        value2 = test_case.values.vector{2,dim2};

        reference = test_case.references.inner.kron{dim1,dim2};

        test_case.evaluate_mtimes("kron",symb1,symb2,value1,value2,reference);
    end
end

methods (Test, ParameterCombination="pairwise", TestTags=["vector" "outer"])
    function test_mtimes_outer(test_case, symb1, symb2, dim1, dim2)
        % Test matrix multiplication on column-by-row vectors.
        value1 = test_case.values.vector{1,dim1};
        value2 = test_case.values.vector{2,dim2}';

        reference = test_case.references.outer.mtimes{dim1,dim2};

        test_case.evaluate_mtimes("mtimes",symb1,symb2,value1,value2,reference);
    end

    function test_kron_outer(test_case, symb1, symb2, dim1, dim2)
        % Test Kronecker product on column-by-row vectors.
        value1 = test_case.values.vector{1,dim1};
        value2 = test_case.values.vector{2,dim2}';

        reference = test_case.references.outer.kron{dim1,dim2};

        test_case.evaluate_mtimes("kron",symb1,symb2,value1,value2,reference);
    end
end

methods (Test, ParameterCombination="pairwise", TestTags="matrix")
    function test_mtimes_matrix(test_case, symb1, symb2, dim1, dim2)
        % Test matrix multiplication on matrices.
        value1 = test_case.values.matrix{1,dim1};
        value2 = test_case.values.matrix{3,dim2};

        reference = test_case.references.matrix.mtimes{dim1,dim2};

        test_case.evaluate_mtimes("mtimes",symb1,symb2,value1,value2,reference);
    end

    function test_kron_matrix(test_case, symb1, symb2, dim1, dim2)
        % Test Kronecker product on matrices.
        value1 = test_case.values.matrix{1,dim1};
        value2 = test_case.values.matrix{3,dim2};

        reference = test_case.references.matrix.kron{dim1,dim2};

        test_case.evaluate_mtimes("kron",symb1,symb2,value1,value2,reference);
    end
end

methods
    function evaluate_mtimes(test_case, op, symb1, symb2, value1, value2, reference)
        % Evaluate matrix operation.
        test_case.assumeTrue(symb1 || symb2, "Operands not symbolic.");
        
        % symbolic polynomials
        [p1,symbol1,argument1] = test_case.get_operand(symb1,value1);
        [p2,symbol2,argument2] = test_case.get_operand(symb2,value2);

        if strcmp(op,"mtimes") && ~any([
                isscalar(value1)
                isscalar(value2)
                size(value1,2) == size(value2,1)
            ])
            % inner dimension mismatch for matrix multiplication
            diagtext = sprintf('Inner dimension mismatch expected: %d vs. %d.',size(value1,2),size(value2,1));
            test_case.verifyError(@() feval(op,p1,p2),?casos.package.core.IncompatibleSizesError,diagtext);
            return
        end

        % else:
        % build symbolic function
        expression = feval(op,p1,p2);
        f = casos.Function('f',[symbol1 symbol2],{expression});

        actual = f(argument1{:},argument2{:});

        % perform assertion
        test_case.verifyEqualPolynomial(actual,reference,"RelTol",1e-15);
    end
end

end
