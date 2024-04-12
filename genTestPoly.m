function  [testValue,refValue] = genTestPoly()


            m    = 3;
            n    = 1;
            deg  = 1:2;
            deg2 = 2;
            
            % generate casos monomial vector
            x_cas        = casos.PS('x',m,n);
            x_monom1_cas = monomials(x_cas,deg);
            x_monom2_cas = monomials(x_cas,deg2);
            
            % generate random coefficients
            coeffs1 = randn(length(x_monom1_cas),1)';
            coeffs2 = randn(length(x_monom2_cas),1)';
            

            x_sopt        = mpvar('x',m,n);
            x_monom1_sopt = monomials(x_sopt,deg);
            x_monom2_sopt = monomials(x_sopt,deg2);
            

            % generate casos polynomials
            poly1_cas = coeffs1 * x_monom1_cas;
            poly2_cas = coeffs2 * x_monom2_cas;
            
            poly1_sopt = coeffs1 * x_monom1_sopt;
            poly2_sopt = coeffs2 * x_monom2_sopt;


            % add test values to testValue
            testValue = [poly1_cas;poly2_cas];

            % add test values to testValue
            refValue = [poly1_sopt;poly2_sopt];

end