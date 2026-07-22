% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

classdef IsEqualPolynomialTo < matlab.unittest.constraints.Constraint
% Constraint to determine if polynomials are equal.

properties (SetAccess=immutable)
    hasSize;
    hasIndeterminates;
    hasCoefficients;
    hasRows;
    hasColumns;
    hasDegrees;

    reference;
end

methods
    function constraint = IsEqualPolynomialTo(reference,varargin)
        % Create new constraint from reference structure.
        assert(isstruct(reference), 'Reference must be a structure.')

        % constraint for size
        constraint.hasSize = matlab.unittest.constraints.HasSize(reference.sz);

        % constraint for indeterminates
        reference_indets = casos.Indeterminates(reference.indets{:});
        constraint.hasIndeterminates = matlab.unittest.constraints.IsEqualTo(reference_indets);
        
        % constraint for nonzero coefficients
        constraint.hasCoefficients = matlab.unittest.constraints.IsEqualTo(reference.coeffs,varargin{:});

        % constraints for tuplet
        constraint.hasRows = matlab.unittest.constraints.IsEqualTo(reference.i);
        constraint.hasColumns = matlab.unittest.constraints.IsEqualTo(reference.j);
        constraint.hasDegrees = matlab.unittest.constraints.IsEqualTo(reference.degrees);

        % store reference polynomial
        sp = casos.Sparsity.tuplet(reference.sz(1),reference.sz(2),reference.i,reference.j,reference_indets,reference.degrees);
        constraint.reference = casos.PD(sp,reference.coeffs);
    end

    function tf = satisfiedBy(constraint,actual)
        % Check if constraint is satisfied by actual polynomial.
        results = satisfiedByConstraints(constraint,actual);
        
        tf = all(results{1,:}); % access booleans from first row of table) 
    end

    function diagnostic = getDiagnosticFor(constraint,actual)
        % Return diagnostic for actual polynomial.
        [results,actual_values] = satisfiedByConstraints(constraint,actual);

        if all(results{1,:})  % access booleans from first row of table
            % constraint passed
            diagnostic = matlab.unittest.diagnostics.ConstraintDiagnostic;
            diagnostic.Description = 'IsEqualPolynomialTo passed.';

        else
            % constraint failed
            diagnostic = matlab.unittest.diagnostics.ConstraintDiagnostic;
            diagnostic.Description = 'IsEqualPolynomialTo failed.';

            if (~results.size)
                % add diagnostic if size constraint failed
                diagnostic_size = getConstraintDiagnosticFor(constraint.hasSize,actual);
                diagnostic_size.Description = 'Size does not match.';
                diagnostic.addCondition(diagnostic_size);
            end

            if (~results.indets)
                % add diagnostic if indeterminates constraint failed
                diagnostic_indets = getConstraintDiagnosticFor(constraint.hasIndeterminates,actual_values.indets);
                diagnostic_indets.Description = 'Indeterminate variables do not match.';
                diagnostic.addCondition(diagnostic_indets);
            end

            if (~results.coeffs)
                % add diagnostic if coefficient constraint failed
                diagnostic_coeffs = getConstraintDiagnosticFor(constraint.hasCoefficients,actual_values.coeffs);
                diagnostic_coeffs.Description = 'Nonzero coefficients do not match.';
                diagnostic.addCondition(diagnostic_coeffs);
            end

            if (~results.rows)
                % add diagnostic if rows constraint failed
                diagnostic_rows = getConstraintDiagnosticFor(constraint.hasRows,actual_values.i);
                diagnostic_rows.Description = 'Row indices do not match.';
                diagnostic.addCondition(diagnostic_rows);
            end

            if (~results.columns)
                % add diagnostic if columns constraint failed
                diagnostic_columns = getConstraintDiagnosticFor(constraint.hasColumns,actual_values.j);
                diagnostic_columns.Description = 'Column indices do not match.';
                diagnostic.addCondition(diagnostic_columns);
            end

            if (~results.degrees)
                % add diagnostic if degrees constraint failed
                diagnostic_degrees = getConstraintDiagnosticFor(constraint.hasDegrees,actual_values.degrees);
                diagnostic_degrees.Description = 'Degrees do not match.';
                diagnostic.addCondition(diagnostic_degrees);
            end
        end

        diagnostic.DisplayDescription = true;
        diagnostic.DisplayConditions = true;
        diagnostic.DisplayActVal = true;
        diagnostic.DisplayExpVal = true;
        diagnostic.ExpValHeader = 'Expected Polynomial:';
        diagnostic.ActValHeader = 'Actual Polynomial:';
        diagnostic.ActVal = actual;
        diagnostic.ExpVal = constraint.reference;
    end
end

methods (Access=private)
    function [results,actuals] = satisfiedByConstraints(constraint,actual)
        % Check if sub-constraints are satisfied by actual polynomial.
        results = table;
        actuals  = struct;

        % only compare nonzero terms
        actual_sparse = sparsify(actual);
        % get indeterminates
        actuals.indets = actual_sparse.indeterminates;
        % get nonzero coefficients
        actuals.coeffs = full(poly2basis(actual_sparse));
        % get tuplet
        [actuals.i,actuals.j,actuals.degrees] = get_tuplet(sparsity(actual_sparse));

        results.size = satisfiedBy(constraint.hasSize,actual);
        results.indets = satisfiedBy(constraint.hasIndeterminates,actuals.indets);
        results.coeffs = satisfiedBy(constraint.hasCoefficients,actuals.coeffs);
        results.rows = satisfiedBy(constraint.hasRows,actuals.i);
        results.columns = satisfiedBy(constraint.hasColumns,actuals.j);
        results.degrees = satisfiedBy(constraint.hasDegrees,actuals.degrees);
    end
end

end
