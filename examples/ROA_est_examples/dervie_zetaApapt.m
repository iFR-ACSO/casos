clear
clc

syms eta eta_max eta_min zeta_max real
syms a b c real


eq1 = zeta_max == a*eta_min^2+b*eta_min+c;
eq2 =        0 == a*eta_max^2+b*eta_max+c;
eq3 =        0 == 2*a*eta_max+b;


sol = solve([eq1;eq2;eq3],[a;b;c]);

a = sol.a;
b = sol.b;
c = sol.c;



matlabFunction(simplify(a*eta^2+b*eta+c))