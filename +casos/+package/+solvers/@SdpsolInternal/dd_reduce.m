function [sdp,args,map,opts] = dd_reduce(obj, sdp, opts, args)
% Reduce a DD cone program to LP 

% check cones
check_cone(obj.get_cones,opts.Kx,'lin');
check_cone(obj.get_cones,opts.Kc,'lin');
check_cone(obj.get_cones,opts.Kx,'dd');
check_cone(obj.get_cones,opts.Kc,'dd');

% get dimensions of cones in the program decision variables
Nlin = get_dimension(obj.get_cones,opts.Kx,'lin'); nlin0  = sum(Nlin);
Ndd  = get_dimension(obj.get_cones,opts.Kx,'dd');  ndd2   = sum(Ndd.^2);
% (unused, but collected for completeness)
Nlor = get_dimension(obj.get_cones,opts.Kx,'lor'); nlor   = sum(Nlor);
Nrot = get_dimension(obj.get_cones,opts.Kx,'rot'); nrot   = sum(Nrot);
Npsd = get_dimension(obj.get_cones,opts.Kx,'psd'); npsd   = sum(Npsd);
Nsdd = get_dimension(obj.get_cones,opts.Kx,'sdd'); nsdd2  = sum(Nsdd.^2);

% get dimensions of cones in the program constraints
Mlin = get_dimension(obj.get_cones,opts.Kc,'lin');
Mdd  = get_dimension(obj.get_cones,opts.Kc,'dd');
% (unused, but collected for completeness)
Mlor = get_dimension(obj.get_cones,opts.Kc,'lor');
Mrot = get_dimension(obj.get_cones,opts.Kc,'rot');
Mpsd = get_dimension(obj.get_cones,opts.Kc,'psd');
Msdd = get_dimension(obj.get_cones,opts.Kc,'sdd');

% temporarily save sdp.x (original)
x_original = sdp.x;
g_original = sdp.g;

% create zero initial 
num_nlin_x = 0;
num_nlin_g = 0;

% verify DD cones in the constraints and create slack DD variables
M_g = cell(length(Mdd),1);
if isfield(opts.Kc,'dd')
    [sdp,args,M_g,num_nlin_g,dd_index_g,opts] = obj.replaceDDcones(sdp, Mdd, Mlin, args, opts, 'g');
end

% verify DD cones in decision variables
M_x = cell(length(Ndd),1);
if isfield(opts.Kx, 'dd')
    [sdp,args,M_x,num_nlin_x,dd_index_x, opts] = obj.replaceDDcones(sdp, Ndd, Mlin, args, opts, 'x');
end

% update decision variables
added_vars = [vertcat(M_g{:}); vertcat(M_x{:})];
n_inserted = length(added_vars);
sdp.x = [sdp.x(1:Nlin); added_vars; sdp.x(Nlin+1:end)];

% add linear cones to constraints
Mlin = Mlin + num_nlin_x + num_nlin_g;

% add the variables M to sdp.x
n_added = ndd2 + n_inserted;
Nlin = Nlin + n_added;

args.dd_lbx = [args.dd_lbx; -inf(n_added,1)];
args.dd_ubx = [args.dd_ubx;  inf(n_added,1)];

% update linear variables and constraints
opts.Kx.lin = Nlin;
opts.Kc.lin = Mlin;

% map from new sdp.x to old (only the non-DD variables)
len_x_orig = length(x_original);
len_x_new  = length(sdp.x);
len_g_orig = length(g_original); 
len_g_new  = length(sdp.g);

% original variables occupy:
%  - first nlin0 entries unchanged
%  - then after insertion, the rest continues
idx_xlin  = 1:nlin0;
idx_xrest = (nlin0 + n_inserted + 1) : len_x_new;

% sanity: lengths must match
assert( length(idx_xlin) + length(idx_xrest) == len_x_orig, ...
    'Indexing mismatch building map.x (lengths do not add up).');

idx_original = [idx_xlin, idx_xrest];

% selection matrix (one 1 per row) <=> map.x = jacobian(x_original, sdp.x)
map.x = sparse(1:len_x_orig, idx_original, 1, len_x_orig, len_x_new);

% map from new sdp.g to old
map.g = [speye(len_g_orig), sparse(len_g_orig, len_g_new-len_g_orig)];

% map from new lam_* to original lam
len_non_dd = nlin0 + nlor + nrot + npsd;

% figure out the columns in sdp.x where those non-DD variables now live
cols_xlin = 1:nlin0;                                            % unchanged
cols_xlor = (nlin0 + n_inserted) + (1:nlor);                    % shifted by n_inserted
cols_xrot = (nlin0 + n_inserted + nlor) + (1:nrot);
cols_xpsd = (nlin0 + n_inserted + nlor + nrot) + (1:npsd);

cols_non_dd = [cols_xlin, cols_xlor, cols_xrot, cols_xpsd];

% sanity checks
assert(all(cols_non_dd >= 1 & cols_non_dd <= len_x_new), ...
    'xtemp: some indices out of bounds in sdp.x.');
assert(length(cols_non_dd) == len_non_dd, ...
    'xtemp: mismatch between expected and actual length.');

% build selection matrix (1 in each row at the column where the variable lives)
% equivalent to xtemp = jacobian([xlin; xlor; xrot; xpsd], sdp.x)
xtemp = sparse(1:len_non_dd, cols_non_dd, 1, len_non_dd, len_x_new);

% DD part uses map.g through equality indices
G = sparse(ndd2, size(map.g,2));
n_eq = dd_index_x.num_eq;              % number of equalities
rows = 1:n_eq;                         % choose first n_eq rows (or any specific ones)
cols = dd_index_x.eq_idx;              % columns from eq_idx

% place ones along diagonal
G(sub2ind(size(G), rows, cols)) = 1;

map.lam_x = [xtemp, sparse(size(xtemp, 1), size(map.g,2));
             zeros(ndd2, size(xtemp,2)), G];

map.lam_a = [sparse(size(map.g,1), size(xtemp,2)), map.g];
end
