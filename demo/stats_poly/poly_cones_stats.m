%% Polynomial optimization.
% study and comparison of SOS, DSOS and SDSOS

% number of indeterminates in p
N_min = 2;
N_max = 20;
N = N_min:N_max;

if ~exist('data_stats.mat', 'file') 
    % pre-allocate data
    time_collector = zeros(length(N), 3);
    ratio_collector = zeros(length(N), 2);
    
    % solve SOS, DSOS and SDSOS programs for different dimensions
    for i = 1:length(N)
        [time_collector(i,:), ratio_collector(i,:)] = main_prog(N(i));
    end

    % save the time data
    save('data_stats.mat', 'time_collector', 'ratio_collector')
else
    load('data_stats.mat')
end

% plot the computation time results
figure(1)
semilogy(N_min:N_max, time_collector(:,1), '-red', 'DisplayName','SOS program')
hold on
semilogy(N_min:N_max, time_collector(:,2), '-.blue', 'DisplayName','DSOS program')
semilogy(N_min:N_max, time_collector(:,3), '--green', 'DisplayName','SDSOS program')
xlabel('n')
ylabel('Time [s]')
grid
xticks(N_min:1:N_max);
legend show

% plot the ratio of the minimum attained (with some normalization with the
% dimension)
figure(2)
plot(N_min:N_max, ratio_collector(:,1)'.^(1./(N_min:N_max)), '-.blue', 'DisplayName','DSOS program')
hold on
plot(N_min:N_max, ratio_collector(:,2)'.^(1./(N_min:N_max)), '--green', 'DisplayName','SDSOS program')
xlabel('n')
ylabel('(sos/relax)^{(1/n)}')
grid
xticks(N_min:1:N_max);
legend show

%% SOS builder and solver 
function [time_collector, ratio_collector] = main_prog(N)
    % indeterminate variable
    x = casos.Indeterminates('x', N);
    
    % set seed value
    % randn('state',0);
    
    % Degree 4 homogeneous
    [coeff, sparse] = poly2basis(casos.PS.sym('a',monomials(x,4:4)));

    % Generate random polynomial
    rancoeff = randn(length(coeff),1);
    p = casos.PS(sparse,rancoeff);
    
    % scalar decision variable
    g = casos.PS.sym('g');
    
    % define SOS problem:
    sos = struct('x',g,'f',-g,'g',p-g*(x'*x)^2);
    
    % solve each program
    [g_sos,   f_sos,   time_sos  ] = sos_prog(sos);
    [g_dsos,  f_dsos,  time_dsos ] = dsos_prog(sos);
    [g_sdsos, f_sdsos, time_sdsos] = sdsos_prog(sos);

    % set outputs
    time_collector = [time_sos, time_dsos, time_sdsos];
    ratio_collector = [f_sos/f_dsos, f_sos/f_sdsos];
end


%% Sum-of-squares cone
function [g_sos, f_sos, elapsed_time] = sos_prog(sos)
    % constraint is scalar SOS cone
    opts = struct('Kc',struct('sos',1), 'Kx', struct('lin',1));
    
    % solve by relaxation to SDP
    S = casos.sossol('S','mosek',sos,opts);
    
    % evaluate
    tic
    sol = S();
    elapsed_time = toc;

    % save solution with SOS cone
    g_sos = full(sol.x);
    f_sos = full(sol.f);
end

%% Diagonally dominant sum-of-squares cone
function [g_dsos, f_dsos, elapsed_time] = dsos_prog(sos)
    % constraint is scalar DSOS cone
    opts = struct('Kc',struct('dsos',1), 'Kx', struct('lin',1));
    
    % solve by relaxation to SDP
    S = casos.sossol('S','mosek',sos,opts);

    % evaluate
    tic
    sol = S();
    elapsed_time = toc;

    % save solution with DSOS cone
    g_dsos = full(sol.x);
    f_dsos = full(sol.f);
end

%% Scaled diagonally dominant sum-of-squares cone
function [g_sdsos, f_sdsos, elapsed_time] = sdsos_prog(sos)
    % constraint is scalar SDSOS cone
    opts = struct('Kc',struct('sdsos',1), 'Kx', struct('lin',1));
    
    % solve by relaxation to SDP
    S = casos.sossol('S','mosek',sos,opts);
    % evaluate
    tic
    sol = S();
    elapsed_time = toc;

    % save solution with SOS cone
    g_sdsos = full(sol.x);
    f_sdsos = full(sol.f);
end