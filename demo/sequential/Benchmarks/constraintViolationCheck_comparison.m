
clear
clc

% Number of random polynomial constraints
N = 10;

% ------------------------------------------------------------------------
% Comment on number of constraints:
% 
% Ten is not unlikely. Think about a problem of a constrained
% region-of-attraction estimation. Assume three controls, which leads already
% to six SOS constraints. One dissipation inequality + positivity
% constraint for lyapunov. If we then have a few additional state
% constraints say two  we have ten in total. 
%
% ------------------------------------------------------------------------


% Maximum degree for random polynomials
maxDeg = 6;


% ------------------------------------------------------------------------
% Comment on maximum degree for random polynomials
% 
% Maximum degree of six is also not unrealistic. Assume a simple
% set-inclusion constraint i.e. s*(V-g) + P is SOS. 
% Assume max degree of s is of degree four and V of degree two 
% we have already a degree six polynomial. 
%
% ------------------------------------------------------------------------

% number of states; problem dependent, but limitied here to five 
% otherwise SOS projection runs in memory problems
nxMax = 12;

% Initialize arrays to store computation times for each method
time_sos                = zeros(1, nxMax-1);   
time_signed_distance    = zeros(1, nxMax-1);  
time_polyopt            =  zeros(1, nxMax-1);
time_sampling_1000      = zeros(1, nxMax-1);   
time_sampling_10000     = zeros(1, nxMax-1);  
time_sampling_100000    = zeros(1, nxMax-1); 

% Loop over nx from 2 to nxMax
for nx = 2:nxMax 
    nx
    % Indeterminate variables
    x = casos.Indeterminates('x', nx);

    % Some random polynomial
    mono = monomials(x, 0:maxDeg);

    % Generate N polynomials with same monomials but random coefficients
    p = casos.PD();
    for k = 1:N
        p = [p; casos.PD(mono, randn(mono.nnz, 1))];
    end

    %% Projection onto SOS cone
    % Gram decision variable
    % s = casos.PS.sym('q', grambasis(p));

    % Projection error
    % e = s - p;

    % Min ||s-p||^2 s.t. s is SOS
    % sos = struct('x', s, 'f', dot(e, e));

    % opts = struct('Kx', struct('sos', N));
    % opts.error_on_fail = 0;

    % Solve by relaxation to SDP
    % S = casos.sossol('S', 'mosek', sos, opts);
    % startSD = tic;
    % sol = S();
    % time_sos(nx-1) = toc(startSD);  % Store time for SOS projection

    %% Define signed distance
    opts = [];
    [~, ~, z] = grambasis(p, ones(length(p), 1));

    % Build unit vectors
    base_s0 = gramunit(z);
    r = casos.PS.sym('r', length(p));
    s0 = casos.PD(base_s0);

    % Min sum(r) s.t. s is SOS
    sos = struct('x', r, 'f', sum(r), 'g', p + r .* s0);
    opts = struct('Kx', struct('lin', length(sos.x)), 'Kc', struct('sos', length(sos.g)));
        
    tic
    % Solve by relaxation to SDP
    S = casos.sossol('S', 'mosek', sos, opts);
    toc
    startSD = tic;
    sol = S();
    time_signed_distance(nx-1) = toc(startSD); 


    %% Define sampling approaches
    % assume a simple box; only needed for computation
    a = -1;
    b = 1;

    % Generate 100000 samples
    x_sample_all = a + (b - a) * rand(nx, 100000);
    x_sample_1000 = num2cell(x_sample_all(:, 1:1000), 2);
    x_sample_10000 = num2cell(x_sample_all(:, 1:10000), 2);
    x_sample_all = num2cell(x_sample_all, 2);

    % Sampling with 1000 samples (first 1000 from the 100000)
    
    startSmp = tic;
    pfun = to_function(p);
    values_1000 = full(pfun(x_sample_1000{:}));
    minVal_1000 = min(min(values_1000));
    time_sampling_1000(nx-1) = toc(startSmp);  % Store time for 1000 samples

    % Sampling with 10000 samples (first 10000 from the 100000)
    
    startSmp = tic;
    values_10000 = full(pfun(x_sample_10000{:}));
    minVal_10000 = min(min(values_10000));
    time_sampling_10000(nx-1) = toc(startSmp);  % Store time for 10000 samples

    % Sampling with 100000 samples
    % startSmp = tic;
    % values_100000 = full(pfun(x_sample_all{:}));
    % minVal_100000 = min(min(values_100000));
    % time_sampling_100000(nx-1) = toc(startSmp);  % Store time for 100000 samples

end

%% Plot the results
figure;
plot(2:nxMax , time_sos, 'r-o', 'LineWidth', 1); hold on;
plot(2:nxMax , time_signed_distance, 'g-*', 'LineWidth', 1);
plot(2:nxMax , time_sampling_1000, 'b-+', 'LineWidth', 1);
plot(2:nxMax , time_sampling_10000, 'm-x', 'LineWidth', 1);
plot(2:nxMax , time_sampling_100000, 'c-s', 'LineWidth', 1);

xlabel('n');
ylabel('Computation Time (seconds)');
legend('SOS Projection', 'Signed Distance', 'Sampling (1000)', 'Sampling (10000)', 'Sampling (100000)','Location','northwest');
% title('Computation Time for Different Methods');
grid on;

% Set the y-axis to logarithmic scale
set(gca, 'YScale', 'log');

% 
% cleanfigure();
% matlab2tikz('compConsVio.tex','width','\figW','height','\figH');

