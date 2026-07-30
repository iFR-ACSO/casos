% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

classdef (TestTags="PS") TestConcatPS < TestSymbolicOperations
% Test concatenation operations on symbolic polynomials.

properties (SetAccess=protected)
    values       % test polynomials
    references   % reference solutions
end

properties (TestParameter)
    nrow = num2cell(1:5);   % numbers of rows to concat
    arg = {1 2};
    symb1 = {false true};
    symb2 = {false true};

    dim1 = num2cell(1:6);
    dim2 = num2cell(1:6);
end

methods (TestClassSetup)
    function initializeTestData(test_case)
        % Initialize test data for unary operations.
        test_case.loadTestData("concat");

        test_case.fatalAssertLength(test_case.nrow,size(test_case.references.scalar.concat,2));
        test_case.fatalAssertLength(test_case.dim1,size(test_case.references.matrix.horzcat,1));
        test_case.fatalAssertLength(test_case.dim2,size(test_case.references.matrix.horzcat,2));
        test_case.fatalAssertLength(test_case.arg,size(test_case.references.scalar.concat,1));
    end
end

methods (Test, ParameterCombination="pairwise", TestTags="scalar")
    function test_concat_scalar(test_case, nrow, arg)
        % Test concatenation to m-by-n matrix.
        ncol = 6-nrow;
        all_values = test_case.values.scalar(arg,1:(nrow*ncol));

        % symbolic polynomial
        [p,symbols,arguments] = cellfun(@(val) test_case.get_operand(true,val), all_values, 'UniformOutput', false);
        symbol = [symbols{:}];
        argument = [arguments{:}];

        args_to_vertcat = mat2cell(p, 1, repmat(nrow,1,ncol));
        args_to_horzcat = cellfun(@(c) vertcat(c{:}), args_to_vertcat, 'UniformOutput', false);

        % build symbolic function
        expression = horzcat(args_to_horzcat{:});
        f = casos.Function('f',symbol,{expression});

        actual = f(argument{:});
        reference = test_case.references.scalar.concat{arg,nrow};

        test_case.verifyEqualPolynomial(actual, reference);
    end
end

methods (Test, ParameterCombination="pairwise", TestTags=["vector" "column"])
    % function test_diagcat_column(test_case, dim1, dim2)
    %     % Test diagonal concatenation of column vectors.
    %     val1 = test_case.values.vector{1,dim1};
    %     val2 = test_case.values.vector{2,dim2};
    % 
    %     actual1 = diagcat(val1,val2);
    %     actual2 = cat(0,val1,val2);
    %     reference = test_case.references.column.diagcat{dim1,dim2};
    % 
    %     test_case.verifyEqualPolynomial(actual1, reference);
    %     test_case.verifyEqualPolynomial(actual2, reference);
    % end

    function test_horzcat_column(test_case, symb1, symb2, dim1, dim2)
        % Test horizontal concatenation of column vectors.
        value1 = test_case.values.vector{1,dim1};
        value2 = test_case.values.vector{2,dim2};

        reference = test_case.references.column.horzcat{dim1,dim2};

        test_case.evaluate_horzcat(symb1,symb2,value1,value2,reference);
    end

    function test_vertcat_column(test_case, symb1, symb2, dim1, dim2)
        % Test vertical concatenation of column vectors.
        value1 = test_case.values.vector{1,dim1};
        value2 = test_case.values.vector{2,dim2};

        reference = test_case.references.column.vertcat{dim1,dim2};

        test_case.evaluate_vertcat(symb1,symb2,value1,value2,reference);
    end
end

methods (Test, ParameterCombination="pairwise", TestTags=["vector" "row"])
    % function test_diagcat_row(test_case, dim1, dim2)
    %     % Test diagonal concatenation of row vectors.
    %     val1 = test_case.values.vector{1,dim1}';
    %     val2 = test_case.values.vector{2,dim2}';
    % 
    %     actual1 = diagcat(val1,val2);
    %     actual2 = cat(0,val1,val2);
    %     reference = test_case.references.row.diagcat{dim1,dim2};
    % 
    %     test_case.verifyEqualPolynomial(actual1, reference);
    %     test_case.verifyEqualPolynomial(actual2, reference);
    % end

    function test_horzcat_row(test_case, symb1, symb2, dim1, dim2)
        % Test horizontal concatenation of row vectors.
        value1 = test_case.values.vector{1,dim1}';
        value2 = test_case.values.vector{2,dim2}';

        reference = test_case.references.row.horzcat{dim1,dim2};

        test_case.evaluate_horzcat(symb1,symb2,value1,value2,reference);
    end

    function test_vertcat_row(test_case, symb1, symb2, dim1, dim2)
        % Test vertical concatenation of row vectors.
        value1 = test_case.values.vector{1,dim1}';
        value2 = test_case.values.vector{2,dim2}';

        reference = test_case.references.row.vertcat{dim1,dim2};
        
        test_case.evaluate_vertcat(symb1,symb2,value1,value2,reference);
    end
end

methods (Test, ParameterCombination="pairwise", TestTags="matrix")
    % function test_diagcat_matrix(test_case, dim1, dim2)
    %     % Test diagonal concatenation of matrices.
    %     val1 = test_case.values.matrix{1,dim1};
    %     val2 = test_case.values.matrix{2,dim2};
    % 
    %     actual1 = diagcat(val1,val2);
    %     actual2 = cat(0,val1,val2);
    %     reference = test_case.references.matrix.diagcat{dim1,dim2};
    % 
    %     test_case.verifyEqualPolynomial(actual1, reference);
    %     test_case.verifyEqualPolynomial(actual2, reference);
    % end

    function test_horzcat_matrix(test_case, symb1, symb2, dim1, dim2)
        % Test horizontal concatenation of matrices.
        value1 = test_case.values.matrix{1,dim1};
        value2 = test_case.values.matrix{2,dim2};

        reference = test_case.references.matrix.horzcat{dim1,dim2};

        test_case.evaluate_horzcat(symb1,symb2,value1,value2,reference);
    end

    function test_vertcat_matrix(test_case, symb1, symb2, dim1, dim2)
        % Test vertical concatenation of matrices.
        value1 = test_case.values.matrix{1,dim1};
        value2 = test_case.values.matrix{2,dim2};

        reference = test_case.references.matrix.vertcat{dim1,dim2};

        test_case.evaluate_vertcat(symb1,symb2,value1,value2,reference);
    end
end

methods
    function evaluate_horzcat(test_case, symb1, symb2, value1, value2, reference)
        % Evaluate horizontal concatenation.
        test_case.assumeTrue(symb1 || symb2, "Operands not symbolic.");
        
        % symbolic polynomials
        [p1,symbol1,argument1] = test_case.get_operand(symb1,value1);
        [p2,symbol2,argument2] = test_case.get_operand(symb2,value2);

        if ~isempty(value1) && ~isempty(value2) && size(value1,1) ~= size(value2,1)
            % first dimension mismatch
            diagtext = sprintf('First dimension mismatch expected: %d vs. %d.',size(value1,1),size(value2,1));
            test_case.verifyError(@() horzcat(p1,p2),?casos.package.core.IncompatibleSizesError,diagtext);
            test_case.verifyError(@() cat(2,p1,p2),?casos.package.core.IncompatibleSizesError,diagtext);
            return
        end

        % else:
        % build symbolic function
        expression1 = horzcat(p1,p2);
        expression2 = cat(2,p1,p2);

        f1 = casos.Function('f',[symbol1 symbol2],{expression1});
        f2 = casos.Function('f',[symbol1 symbol2],{expression2});

        actual1 = f1(argument1{:},argument2{:});
        actual2 = f2(argument1{:},argument2{:});

        test_case.verifyEqualPolynomial(actual1, reference);
        test_case.verifyEqualPolynomial(actual2, reference);
    end

    function evaluate_vertcat(test_case, symb1, symb2, value1, value2, reference)
        % Evaluate vertical concatenation.
        test_case.assumeTrue(symb1 || symb2, "Operands not symbolic.");
        
        % symbolic polynomials
        [p1,symbol1,argument1] = test_case.get_operand(symb1,value1);
        [p2,symbol2,argument2] = test_case.get_operand(symb2,value2);

        if ~isempty(value1) && ~isempty(value2) && size(value1,2) ~= size(value2,2)
            % second dimension mismatch
            diagtext = sprintf('Second dimension mismatch expected: %d vs. %d.',size(value1,2),size(value2,2));
            test_case.verifyError(@() vertcat(p1,p2),?casos.package.core.IncompatibleSizesError,diagtext);
            test_case.verifyError(@() cat(1,p1,p2),?casos.package.core.IncompatibleSizesError,diagtext);
            return
        end

        % else
        % build symbolic function
        expression1 = vertcat(p1,p2);
        expression2 = cat(1,p1,p2);

        f1 = casos.Function('f',[symbol1 symbol2],{expression1});
        f2 = casos.Function('f',[symbol1 symbol2],{expression2});

        actual1 = f1(argument1{:},argument2{:});
        actual2 = f2(argument1{:},argument2{:});

        test_case.verifyEqualPolynomial(actual1, reference);
        test_case.verifyEqualPolynomial(actual2, reference);
    end
end

end
