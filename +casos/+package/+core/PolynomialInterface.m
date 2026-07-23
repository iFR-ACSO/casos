% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

classdef (Abstract) PolynomialInterface < casos.package.core.Printable
% Base class for polynomial-like objects.

methods (Abstract)
    tf = is_wellposed(obj);
end

methods
    %% Concatenation interface
    function p = horzcat(varargin)
        % Horizontal concatenation.
        p = cat(2,varargin{:});
    end

    function p = vertcat(varargin)
        % Vertical concatenation.
        p = cat(1,varargin{:});
    end

    function p = blkdiag(varargin)
        % Block diagonal concatenation.
        p = cat(0,varargin{:});
    end

    function p = cat(dim,varargin)
        % Generic concatenation.
        switch (nargin-1)
            case 0, p = [];
            case 1, p = varargin{:};
            case 2, error('Notify the developers.')
            otherwise
                % more than two inputs
                N = (nargin-1)/2;
                % recursion
                p1 = cat(dim,varargin{1:floor(N)});
                p2 = cat(dim,varargin{floor(N)+1:end}); 
                % concatenate
                p = cat(dim,p1,p2);
        end
    end
end

methods (Access=protected)
    %% Utilities
    function [tf,sz] = check_sz_assign(a,b)
        % Check if sizes are compatible for matrix assignment.
        sza = size(a);
        szb = size(b);

        % dimensions are compatible if equal or right side is row/column
        tf = all(sza == szb | szb == 1);
        % size of left side does not change
        sz = sza;
    end

    function [tf,sz] = check_sz_basic(a,b)
        % Check if sizes are compatible for basic (array) operations.
        sza = size(a);
        szb = size(b);
        
        % compare dimensions
        I = (sza == szb);
        % find zero dimension
        I0 = (sza == 0) | (szb == 0);
        % find one dimension
        I1 = (sza == 1) | (szb == 1);
        
        % dimensions are compatible if equal or one summand is row/column
        tf = all(I | I1);

        % dimensions of result
        sz = max(sza,szb);
        % operation with empty operands is empty
        sz(I0) = 0;
    end

    function [tf,sz] = check_sz_concat(dim,a,b)
        % Check if sizes are compatible for concatenation.
        sza = size(a);
        szb = size(b);

        % find zero dimension
        I0 = (sza == 0) | (szb == 0);

        % dimensions are consistent if size(a,~dim) == size(b,~dim)
        % unless either operand is empty or diagonal concatenation
        tf = (dim == 0) || any(I0) || (sza(3-dim) == (szb(3-dim)));

        if (dim == 0) || (sza(3-dim) == 0 && szb(3-dim) == 0)
            % concatenate both (non-empty) dimensions
            sz = sza + szb;
        elseif any(sza == 0) && any(szb == 0)
            % concatenation of empty polynomials
            sz(dim) = min(sza(dim), szb(dim));
            sz(3-dim) = max(sza(3-dim), szb(3-dim));
        elseif any(sza == 0)
            % concatenation of b with empty polynomial
            sz = szb;
        elseif any(szb == 0)
            % concatenation of a with empty polynomial
            sz = sza;
        else
            % concatenate along dimension
            sz(dim) = sza(dim) + szb(dim);
            sz(3-dim) = sza(3-dim);
        end
    end

    function [tf,sz] = check_sz_equal(a,b)
        % Check if two objects have same size.
        sz = size(a);
        tf = isequal(sz,size(b));
    end

    function [tf,sz] = check_sz_matrix(a,b)
        % Check if sizes are compatible for matrix multiplication.
        tf = size(a,2) == size(b,1);
        sz = [size(a,1) size(b,2)];
    end

    function [tf,sz] = check_sz_square(a)
        % Check if matrix is square.
        sz = size(a);
        tf = isequal(sz(1),sz(2));
    end
end

end
