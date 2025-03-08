
nx = 1;
x  = casos.Indeterminates('x',nx);
y  = casos.Indeterminates('y',nx);
z  = casos.Indeterminates('z',nx);


p = 12+y^2-2*x^3*y+2*y*z^2+x^6-2*x^3*z^2+z^4+x^2*y^2



[~,~,z] = grambasis(p, ones(length(p),1 ) );

% build unit vectors
base_s0 = gramunit(z);

r  = casos.PS.sym('r',length(p));
s0 = casos.PD(base_s0)