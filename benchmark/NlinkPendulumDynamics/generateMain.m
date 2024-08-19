% Test and plot the dynamics of a N-link inverted pendulum

% Number of links
for n = 2:10

    % degree of the polynomial approximation
    deg = 2;
    
    pendulum_dyn_ctrl = ['pendulum_dyn_ctrl_n' num2str(n)];
    pendulum_dyn_poly = ['pendulum_dyn_poly_n' num2str(n) '_d' num2str(deg)];
    
    % generate 
    % 1. full nlink pendulum dynamics
    % 2. get the linerize model at the upright position
    % 3. get polynomial model of degree 2
    data_filename = ['data_n' int2str(n) '.mat'];
    if exist(data_filename, 'file')~=2 || ...
            exist(pendulum_dyn_poly, 'file')~=2
        [K,S] = generate_nlink_pendulum(n,deg);
        save(data_filename, "S","K")
    else
        % load(data_filename)
    end

    disp(['Generated ' num2str(n) ' -link pendulum dynamics and corresponding local control law.'])

end