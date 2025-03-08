% Generate N-link pendulum dynamics


% degree of the polynomial approximation
deg = 2;

% if only the polynomial dynamics are desired, use only_poly=1 
only_poly = 0;

for n = 2:10
pendulum_dyn_ctrl = ['pendulum_dyn_ctrl_n' num2str(n)];
pendulum_dyn_poly = ['pendulum_dyn_poly_n' num2str(n) '_d' num2str(deg)];

% generate 
% 1. full nlink pendulum dynamics
% 2. get the linerize model at the upright position
% 3. get polynomial model of degree 2
data_filename = ['data_n' int2str(n) '.mat'];

% only generate, if it does not exist yet
% if exist(data_filename, 'file')~=2 || ...
        % exist(pendulum_dyn_poly, 'file')~=2
    startTime = tic;
    [K,P,A,B] = generate_nlink_pendulum(n,deg, only_poly);
    save(data_filename, "K","P","A","B")
% else
%     load(data_filename)
% end

disp(['Generated ' num2str(n) '-linked pendulum. Took: ' num2str(toc(startTime)) 's'])

end