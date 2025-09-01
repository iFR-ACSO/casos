function [sdp,args,map,opts] = sdd_reduce(obj, sdp, opts, args)
% Reduce a SDD cone program to SOCP 

% check cones
check_cone(obj.get_cones,opts.Kx,'lin');
check_cone(obj.get_cones,opts.Kc,'lin');
check_cone(obj.get_cones,opts.Kx,'sdd');
check_cone(obj.get_cones,opts.Kc,'sdd');

% get dimensions of cones in the program decision variables
Nlin = get_dimension(obj.get_cones,opts.Kx,'lin'); nlin0  = sum(Nlin);
Nsdd = get_dimension(obj.get_cones,opts.Kx,'sdd'); nsdd2  = sum(Nsdd.^2);
Nlor = get_dimension(obj.get_cones,opts.Kx,'lor'); nlor   = sum(Nlor);
Nrot = get_dimension(obj.get_cones,opts.Kx,'rot'); nrot   = sum(Nrot);
Npsd = get_dimension(obj.get_cones,opts.Kx,'psd'); npsd   = sum(Npsd);
% (unused, but collected for completeness)
Ndd  = get_dimension(obj.get_cones,opts.Kx,'dd');  ndd2   = sum(Ndd.^2);

% get dimensions of cones in the program constraints
Mlin = get_dimension(obj.get_cones,opts.Kc,'lin');
Msdd = get_dimension(obj.get_cones,opts.Kc,'sdd');
Mpsd = get_dimension(obj.get_cones,opts.Kc,'psd');
Mdd  = get_dimension(obj.get_cones,opts.Kc,'dd');
% (unused, but collected for completeness)
Mlor = get_dimension(obj.get_cones,opts.Kc,'lor');
Mrot = get_dimension(obj.get_cones,opts.Kc,'rot');

% save sizes of sdp.x and sdp.g (original)
len_x_orig = length(sdp.x);     % original decision variable length
len_g_orig = length(sdp.g);     % original constraint length

% initialize SDD index structs in case the fields don't exist
sdd_index_g.num_eq = 0;         % in constraints
sdd_index_x.num_eq = 0;         % in decision variables

% verify SDD cones in the constraints and create slack SDD variables
M_g = cell(length(Mdd),1);
if isfield(opts.Kc,'sdd')
    [sdp,args,M_g,~,sdd_index_g,opts] = obj.replaceSDDcones(sdp, Msdd, Mlin, args, opts, 'g');
end

% verify SDD cones in decision variables
M_x = cell(length(Ndd),1);
if isfield(opts.Kx, 'sdd')
    [sdp,args,M_x,~,sdd_index_x, opts] = obj.replaceSDDcones(sdp, Nsdd, Mlin, args, opts, 'x');
end

% update decision variables
added_vars = [vertcat(M_g{:}); vertcat(M_x{:})];
n_inserted = length(added_vars);
sdp.x = [sdp.x(1:Nlin); added_vars; sdp.x(Nlin+1:end)];

% add the variables M to sdp.x
Nlin = Nlin + nsdd2+n_inserted;
Mlin = Mlin + sdd_index_x.num_eq+sdd_index_g.num_eq; 
Mpsd = [Mpsd, repmat(2, 1, n_inserted/3)];

% update lbx and ubx 
args.dd_lbx = [args.dd_lbx; -inf(nsdd2 + n_inserted,1)];
args.dd_ubx = [args.dd_ubx;  inf(nsdd2 + n_inserted,1)];

% update cones
opts.Kx.lin = Nlin;
opts.Kc.lin = Mlin;
opts.Kc.psd = Mpsd;

% get sizes of NEW sdp.x and sdp.g
len_x_new  = length(sdp.x);
len_g_new  = length(sdp.g);

% original variables occupy:
%  - first nlin0 entries unchanged
%  - then after insertion, the rest continues
idx_xlin  = 1:nlin0;
idx_xrest = (nlin0 + n_inserted + 1) : len_x_new;

% sanity check
assert(length(idx_xlin) + length(idx_xrest) == len_x_orig, ...
    'Indexing mismatch building map.x (lengths do not add up).');

idx_original = [idx_xlin, idx_xrest];

% build sparse selection matrix: one 1 per row
% equivalent to map.x = jacobian(x_original, sdp.x);
map.x = sparse(1:len_x_orig, idx_original, 1, len_x_orig, len_x_new);

% map from new sdp.g to old
map.g = [speye(len_g_orig), sparse(len_g_orig, len_g_new-len_g_orig)];

% map.lam_x for non-SDD variables
len_non_sdd = nlin0 + nlor + nrot + npsd;

% figure out the columns in sdp.x where non-SDD variables now live
cols_xlin = 1:nlin0;                                          % unchanged
cols_xlor = (nlin0 + n_inserted) + (1:nlor);                  % shifted by inserted vars
cols_xrot = (nlin0 + n_inserted + nlor) + (1:nrot);
cols_xpsd = (nlin0 + n_inserted + nlor + nrot) + (1:npsd);

cols_non_sdd = [cols_xlin, cols_xlor, cols_xrot, cols_xpsd];

% sanity checks
assert(all(cols_non_sdd >= 1 & cols_non_sdd <= len_x_new), ...
    'xtemp: some indices out of bounds in sdp.x.');
assert(length(cols_non_sdd) == len_non_sdd, ...
    'xtemp: mismatch between expected and actual length.');

% selection matrix for non-SDD variables
% equivalent to xtemp = jacobian([xlin; xlor; xrot; xpsd; xdd], sdp.x);
map_non_sdd = sparse(1:len_non_sdd, cols_non_sdd, 1, len_non_sdd, len_x_new);

% map.lam_x for SDD variables using equality indices
map_sdd_eqs = sparse(nsdd2, len_g_new);
n_eq = sdd_index_x.num_eq;              % number of equalities
rows = 1:n_eq;                          % choose first n_eq rows (or any specific ones)
cols = sdd_index_x.eq_idx;              % columns from eq_idx

% place ones along diagonal
map_sdd_eqs(sub2ind([nsdd2, len_g_new], rows, cols)) = 1;

% combine non-SDD and SDD selection matrices
map.lam = [map_non_sdd, sparse(len_non_sdd, len_g_new);
             sparse(nsdd2, len_x_new), map_sdd_eqs];

map.lam = [map.lam;
           sparse(len_g_orig, len_x_new), map.g];

% map lam_a (dual variables for g)
% map.lam_a = [sparse(len_g_orig, len_x_new), map.g];

end




