% Test Newton polytope monomial basis reduction 

x = casos.PS('x', 3, 1);
f = casos.PS.sym('f');

g = [12+x(2)^2-2*x(1)^3*x(2)+2*x(2)*x(3)^2+x(1)^6-2*x(1)^3*x(3)^2+x(3)^4+x(1)^2*x(2)^2 + f;
    (x(1)-1)^2+x(1)^4+f;
    3*x(1)^4-2*x(1)^2*x(2)+7*x(1)^2-4*x(1)*x(2)+4*x(2)^2+f;
    x(1)^4 + x(2)^4 - 2*x(1)^2*x(2)^2];

sos = struct('x', f, 'g', g, 'f', f);
% states + constraint are SOS cones
opts.Kx.l = 1;
opts.Kx.s = 0;
opts.Kc.s = 4;

% ignore infeasibility
opts.error_on_fail = false;

% solve by relaxation to SDP
S = casos.sossol('S','sedumi',sos,opts);

N_times = 1000;
idx = 1:N_times;

opts.newton = 0; tic;
arrayfun(@(i) S(), idx, 'UniformOutput', false);
fprintf('newton = %d: %ds \n', opts.newton ,toc)

opts.newton = 1; tic;
arrayfun(@(i) S(), idx, 'UniformOutput', false);
fprintf('newton = %d: %ds \n', opts.newton ,toc)



