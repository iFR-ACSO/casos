% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

classdef (TestTags="PS") TestSubsPS < TestSymbolicOperations
% Test substitution in symbolic polynomials.

properties (Constant)
    nnz_threshold = 100;
end

properties (SetAccess=protected)
    values       % test polynomials
    references   % reference solutions
end

properties (TestParameter)
    ivar = num2cell(1:4);
    symb1 = {false true};
    symb2 = {false true};

    dim =  num2cell(1:6);
    arg = num2cell(1:10);
end

methods (TestClassSetup)
    function initializeTestData(test_case)
        % Initialize test data for substitution operation.
        test_case.loadTestData("subs");

        test_case.fatalAssertLength(test_case.ivar,size(test_case.references.scalar.single,2));
        test_case.fatalAssertLength(test_case.dim,size(test_case.references.matrix.single,1));
        test_case.fatalAssertLength(test_case.arg,size(test_case.references.scalar.single,1));
    end
end

methods (Test, ParameterCombination="pairwise", TestTags="scalar")
    function test_subs_single_scalar(test_case, symb1, symb2, ivar, arg)
        % Test substitution of a single variable.
        value1 = test_case.values.scalar{1,arg};
        value2 = test_case.values.scalar{2,arg};

        % reduce order of polynomials
        value1 = project(value1,restrict_terms(sparsity(value1),0:4));
        value2 = project(value2,restrict_terms(sparsity(value2),0:4));

        reference = test_case.references.scalar.single{arg,ivar};

        test_case.evaluate_subs_single(symb1,symb2,ivar,value1,value2,reference);
    end

    function test_subs_multiple_scalar(test_case, symb1, symb2, ivar, arg)
        % Test substitution of multiple variables.
        value1 = test_case.values.scalar{2,arg};
        value2 = test_case.values.vector{2,value1.nvars+1};

        % reduce order of polynomials
        value1 = project(value1,restrict_terms(sparsity(value1),0:4));
        value2 = project(value2,restrict_terms(sparsity(value2),0:4));

        reference = test_case.references.scalar.multiple{arg,ivar};

        test_case.evaluate_subs_multiple(symb1,symb2,ivar,value1,value2,reference);
    end
end

methods (Test, ParameterCombination="pairwise", TestTags=["vector" "column"])
    function test_subs_single_column(test_case, symb1, symb2, ivar, dim)
        % Test substitution of a single variable in a column vector.
        value1 = test_case.values.vector{1,dim};
        value2 = test_case.values.scalar{1,dim};

        % reduce order of polynomials
        value1 = project(value1,restrict_terms(sparsity(value1),0:4));
        value2 = project(value2,restrict_terms(sparsity(value2),0:4));

        reference = test_case.references.column.single{dim,ivar};

        test_case.evaluate_subs_single(symb1,symb2,ivar,value1,value2,reference);
    end

    function test_subs_multiple_column(test_case, symb1, symb2, ivar, dim)
        % Test substitution of multiple variables in a column vector.
        value1 = test_case.values.vector{2,dim};
        value2 = test_case.values.vector{1,value1.nvars+1};

        % reduce order of polynomials
        value1 = project(value1,restrict_terms(sparsity(value1),0:4));
        value2 = project(value2,restrict_terms(sparsity(value2),0:4));

        reference = test_case.references.column.multiple{dim,ivar};

        test_case.evaluate_subs_multiple(symb1,symb2,ivar,value1,value2,reference);
    end
end

methods (Test, ParameterCombination="pairwise", TestTags=["vector" "row"])
    function test_subs_single_row(test_case, symb1, symb2, ivar, dim)
        % Test substitution of a single variable in a row vector.
        value1 = test_case.values.vector{2,dim}';
        value2 = test_case.values.scalar{2,dim};

        % reduce order of polynomials
        value1 = project(value1,restrict_terms(sparsity(value1),0:4));
        value2 = project(value2,restrict_terms(sparsity(value2),0:4));

        reference = test_case.references.row.single{dim,ivar};

        test_case.evaluate_subs_single(symb1,symb2,ivar,value1,value2,reference);
    end

    function test_subs_multiple_row(test_case, symb1, symb2, ivar, dim)
        % Test substitution of multiple variables in a row vector.
        value1 = test_case.values.vector{1,dim}';
        value2 = test_case.values.vector{2,value1.nvars+1};

        % reduce order of polynomials
        value1 = project(value1,restrict_terms(sparsity(value1),0:4));
        value2 = project(value2,restrict_terms(sparsity(value2),0:4));

        reference = test_case.references.row.multiple{dim,ivar};

        test_case.evaluate_subs_multiple(symb1,symb2,ivar,value1,value2,reference);
    end
end

methods (Test, ParameterCombination="pairwise", TestTags="matrix")
    function test_subs_single_matrix(test_case, symb1, symb2, ivar, dim)
        % Test substitution of a single variable in a matrix.
        value1 = test_case.values.matrix{2,dim};
        value2 = test_case.values.scalar{1,dim};

        % reduce order of polynomials
        value1 = project(value1,restrict_terms(sparsity(value1),0:4));
        value2 = project(value2,restrict_terms(sparsity(value2),0:4));

        reference = test_case.references.matrix.single{dim,ivar};

        test_case.evaluate_subs_single(symb1,symb2,ivar,value1,value2,reference);
    end

    function test_subs_multiple_matrix(test_case, symb1, symb2, ivar, dim)
        % Test substitution of multiple variables in a column vector.
        value1 = test_case.values.matrix{1,dim};
        value2 = test_case.values.vector{1,value1.nvars+1};

        % reduce order of polynomials
        value1 = project(value1,restrict_terms(sparsity(value1),0:4));
        value2 = project(value2,restrict_terms(sparsity(value2),0:4));

        reference = test_case.references.matrix.multiple{dim,ivar};

        test_case.evaluate_subs_multiple(symb1,symb2,ivar,value1,value2,reference);
    end
end

methods
    function evaluate_subs_single(test_case, symb1, symb2, ivar, value1, value2, reference)
        % Evaluate substitution of single variable.
        test_case.assumeTrue(symb1 || symb2, "Operands not symbolic.");
        
        test_case.assumeLessThan(nnz(value1),test_case.nnz_threshold,'Number of nonzero coefficients exceeds threshold.')
        test_case.assumeLessThan(nnz(value2),test_case.nnz_threshold,'Number of nonzero coefficients exceeds threshold.')

        vars = value1.indeterminates;
        if ivar > length(vars)
            % variable not in polynomial
            var = casos.Indeterminates('y');
        else
            var = vars(ivar);
        end

        % symbolic polynomial
        [p1,symbol1,argument1] = test_case.get_operand(symb1,value1);
        [p2,symbol2,argument2] = test_case.get_operand(symb2,value2);

        % build symbolic function
        expression = subs(p1,var,p2);
        f = casos.Function('f',[symbol1 symbol2],{expression});

        actual = f(argument1{:},argument2{:});

        % perform assertion
        test_case.verifyEqualPolynomial(actual,reference,"RelTol",1e-12);
    end

    function evaluate_subs_multiple(test_case, symb1, symb2, ivar, value1, value2, reference)
        % Evaluate substitution of multiple variables.
        test_case.assumeTrue(symb1 || symb2, "Operands not symbolic.");
        
        test_case.assumeLessThan(nnz(value1),test_case.nnz_threshold,'Number of nonzero coefficients exceeds threshold.')
        test_case.assumeLessThan(nnz(value2),test_case.nnz_threshold,'Number of nonzero coefficients exceeds threshold.')

        vars = value1.indeterminates;
        if ivar <= length(vars)
            % replace variable
            vars(ivar) = casos.Indeterminates('y');
        end
        
        % symbolic polynomial
        [p1,symbol1,argument1] = test_case.get_operand(symb1,value1);
        [p2,symbol2,argument2] = test_case.get_operand(symb2,value2);

        % build symbolic function
        expression = subs(p1,vars,p2);
        f = casos.Function('f',[symbol1 symbol2],{expression});

        actual = f(argument1{:},argument2{:});

        % perform assertion
        test_case.verifyEqualPolynomial(actual,reference,"RelTol",1e-12);
    end
end

end
