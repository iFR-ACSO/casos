% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

classdef (TestTags="PD") TestConcatPD < TestPolynomialOperations
% Test concatenation operations on constant polynomials.

properties (SetAccess=protected)
    values       % test polynomials
    references   % reference solutions
end

properties (TestParameter)
    nrow = num2cell(1:5);   % numbers of rows to concat
    arg = {1 2};

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
        value = test_case.values.scalar(arg,1:(nrow*ncol));

        args_to_vertcat = mat2cell(value, 1, repmat(nrow,1,ncol));
        args_to_horzcat = cellfun(@(c) vertcat(c{:}), args_to_vertcat, 'UniformOutput', false);

        actual = horzcat(args_to_horzcat{:});
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

    function test_horzcat_column(test_case, dim1, dim2)
        % Test horizontal concatenation of column vectors.
        value1 = test_case.values.vector{1,dim1};
        value2 = test_case.values.vector{2,dim2};

        reference = test_case.references.column.horzcat{dim1,dim2};

        test_case.evaluate_horzcat(value1,value2,reference);
    end

    function test_vertcat_column(test_case, dim1, dim2)
        % Test vertical concatenation of column vectors.
        value1 = test_case.values.vector{1,dim1};
        value2 = test_case.values.vector{2,dim2};

        reference = test_case.references.column.vertcat{dim1,dim2};

        test_case.evaluate_vertcat(value1,value2,reference);
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

    function test_horzcat_row(test_case, dim1, dim2)
        % Test horizontal concatenation of row vectors.
        value1 = test_case.values.vector{1,dim1}';
        value2 = test_case.values.vector{2,dim2}';

        reference = test_case.references.row.horzcat{dim1,dim2};

        test_case.evaluate_horzcat(value1,value2,reference);
    end

    function test_vertcat_row(test_case, dim1, dim2)
        % Test vertical concatenation of row vectors.
        value1 = test_case.values.vector{1,dim1}';
        value2 = test_case.values.vector{2,dim2}';

        reference = test_case.references.row.vertcat{dim1,dim2};
        
        test_case.evaluate_vertcat(value1,value2,reference);
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

    function test_horzcat_matrix(test_case, dim1, dim2)
        % Test horizontal concatenation of matrices.
        value1 = test_case.values.matrix{1,dim1};
        value2 = test_case.values.matrix{2,dim2};

        reference = test_case.references.matrix.horzcat{dim1,dim2};

        test_case.evaluate_horzcat(value1,value2,reference);
    end

    function test_vertcat_matrix(test_case, dim1, dim2)
        % Test vertical concatenation of matrices.
        value1 = test_case.values.matrix{1,dim1};
        value2 = test_case.values.matrix{2,dim2};

        reference = test_case.references.matrix.vertcat{dim1,dim2};

        test_case.evaluate_vertcat(value1,value2,reference);
    end
end

methods
    function evaluate_horzcat(test_case, value1, value2, reference)
        % Evaluate horizontal concatenation.
        if ~isempty(value1) && ~isempty(value2) && size(value1,1) ~= size(value2,1)
            % first dimension mismatch
            diagtext = sprintf('First dimension mismatch expected: %d vs. %d.',size(value1,1),size(value2,1));
            test_case.verifyError(@() horzcat(value1,value2),?casos.package.core.IncompatibleSizesError,diagtext);
            test_case.verifyError(@() cat(2,value1,value2),?casos.package.core.IncompatibleSizesError,diagtext);
            return
        end

        % else
        actual1 = horzcat(value1,value2);
        actual2 = cat(2,value1,value2);

        test_case.verifyEqualPolynomial(actual1, reference);
        test_case.verifyEqualPolynomial(actual2, reference);
    end

    function evaluate_vertcat(test_case, value1, value2, reference)
        % Evaluate vertical concatenation.
        if ~isempty(value1) && ~isempty(value2) && size(value1,2) ~= size(value2,2)
            % second dimension mismatch
            diagtext = sprintf('Second dimension mismatch expected: %d vs. %d.',size(value1,2),size(value2,2));
            test_case.verifyError(@() vertcat(value1,value2),?casos.package.core.IncompatibleSizesError,diagtext);
            test_case.verifyError(@() cat(1,value1,value2),?casos.package.core.IncompatibleSizesError,diagtext);
            return
        end

        % else
        actual1 = vertcat(value1,value2);
        actual2 = cat(1,value1,value2);

        test_case.verifyEqualPolynomial(actual1, reference);
        test_case.verifyEqualPolynomial(actual2, reference);
    end
end

end
