function  [testValue,refValue] = genTestPoly()

            lower = 1;
            upper = 4;

            m    = ceil(lower + (upper-lower)*rand());  % number of indeterminates

            deg  = 0:ceil(lower + (upper-lower)*rand()); % degree of first polynomial
            deg2 = 1:ceil(lower + (upper-lower)*rand());   % degree of second polynomial
            
            % generate casos monomial vector
            x_cas        = casos.PS('x',m,1);
            x_monom1_cas = monomials(x_cas,deg);
            x_monom2_cas = monomials(x_cas,deg2);
            
            % generate random coefficients
            coeffs1 = rand(x_monom1_cas.nnz,1)';
            coeffs2 = rand(x_monom2_cas.nnz,1)';
            

            % generate casos polynomials
            poly1_cas = casos.PS(x_monom1_cas,coeffs1);
            poly2_cas = casos.PS(x_monom2_cas,coeffs2);

            % generate reference polynomials using multipoly
            x_sopt        = mpvar('x',m,1);
            x_monom1_sopt = monomials(x_sopt,deg);
            x_monom2_sopt = monomials(x_sopt,deg2);
            
          
            poly1_sopt = coeffs1 * x_monom1_sopt;
            poly2_sopt = coeffs2 * x_monom2_sopt;

            % add test values to testValue
            testValue = [poly1_cas;poly2_cas];

            % add test values to testValue
            refValue = [poly1_sopt;poly2_sopt];

end