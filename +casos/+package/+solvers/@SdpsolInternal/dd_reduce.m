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
Nlor = get_dimension(obj.get_cones,opts.Kx,'lor'); nlor   = sum(Nlor);
Nrot = get_dimension(obj.get_cones,opts.Kx,'rot'); nrot   = sum(Nrot);
Npsd = get_dimension(obj.get_cones,opts.Kx,'psd'); npsd   = sum(Npsd.^2);

% get dimensions of cones in the program constraints
Mlin = get_dimension(obj.get_cones,opts.Kc,'lin'); mlin0 = sum(Mlin);
Mdd  = get_dimension(obj.get_cones,opts.Kc,'dd');  mdd2  = sum(Mdd.^2);
Mlor = get_dimension(obj.get_cones,opts.Kc,'lor'); mlor  = sum(Mlor);
Mrot = get_dimension(obj.get_cones,opts.Kc,'rot'); mrot  = sum(Mrot);
Mpsd = get_dimension(obj.get_cones,opts.Kc,'psd'); mpsd  = sum(Mpsd.^2);

% save sizes of sdp.x and sdp.g (original)
len_x_orig = size(sdp.x, 1);     % original decision variable length
len_g_orig = size(sdp.g, 1);     % original constraint length

% create zero initial 
num_nlin_x = 0;     % in decision variables
num_nlin_g = 0;     % in constraints

% verify DD cones in the constraints and create slack DD variables
M_g = cell(length(Mdd),1);
if mdd2 > 0
    [sdp,args,M_g,num_nlin_g,dd_index_g,opts] = replaceDDcones(sdp, Mdd, Mlin, args, opts, 'g');
else
    dd_index_g.num_eq = 0;
    dd_index_g.eq_idx = [];
end

% verify DD cones in decision variables
M_x = cell(length(Ndd),1);
if ndd2 > 0
    [sdp,args,M_x,num_nlin_x,dd_index_x, opts] = replaceDDcones(sdp, Ndd, Mlin, args, opts, 'x');
else
    dd_index_x.num_eq = 0;
    dd_index_x.eq_idx = [];
end

% update decision variables
added_vars = [vertcat(M_g{:}); vertcat(M_x{:})];
n_inserted = length(added_vars);

% move the DD variables to the lin position
ddvar  = sdp.x(nlin0+nlor+nrot+npsd+1:end);
psdvar = sdp.x(nlin0+nlor+nrot+1:nlin0+nlor+nrot+npsd);
sdp.x  = [sdp.x(1:Nlin); added_vars; ddvar; psdvar]; 

% add linear cones to constraints
m_inserted = num_nlin_x + num_nlin_g;
Mlin = Mlin + m_inserted;

% add the variables M to sdp.x
n_added = ndd2 + n_inserted;
Nlin = Nlin + n_added;

args.dd_lbx = [args.dd_lbx; -inf(n_added,1)];
args.dd_ubx = [args.dd_ubx;  inf(n_added,1)];

% update cones
opts.Kx.lin = Nlin;
opts.Kc.lin = Mlin;

% get sizes of NEW sdp.x and sdp.g
len_x_new  = length(sdp.x);
len_g_new  = length(sdp.g);

% original variables occupy:
%  - first nlin0 entries unchanged
%  - then after insertion, the rest continues
idx_xlin  = 1:nlin0;
idx_xpsd  = (nlin0+1):(nlin0+npsd);
idx_xdd   = (nlin0+npsd+1):(nlin0+npsd+ndd2);

idx_xrest = (nlin0 + n_inserted + 1) : len_x_new;

% sanity: lengths must match
assert( length(idx_xlin) + length(idx_xrest) == len_x_orig, ...
    'Indexing mismatch building map.x (lengths do not add up).');

%idx_original = [idx_xlin, n_inserted+idx_xdd, n_inserted+idx_xpsd];
idx_original = [idx_xlin, ndd2+n_inserted+idx_xpsd, n_inserted+idx_xdd-npsd];


% selection matrix (one 1 per row) <=> map.x = jacobian(x_original, sdp.x)
map.x = sparse(1:len_x_orig, idx_original, 1, len_x_orig, len_x_new);

% map from new sdp.g to old
map.g = [speye(len_g_orig), sparse(len_g_orig, len_g_new-len_g_orig)];

% map from new lam_* to original lam
len_non_ndd = nlin0 + nlor + nrot + npsd;
len_non_mdd = mlin0 + mlor + mrot + mpsd;

% figure out the columns in sdp.x where those non-DD variables now live
cols_xlin = 1:nlin0;                                            % unchanged
cols_xlor = (nlin0 + n_inserted) + (1:nlor);                    % shifted by n_inserted
cols_xrot = (nlin0 + n_inserted + nlor) + (1:nrot);
cols_xpsd = (nlin0 + n_inserted + nlor + nrot) + (1:npsd);

cols_non_dd = [cols_xlin, cols_xlor, cols_xrot, cols_xpsd];

% sanity checks
assert(all(cols_non_dd >= 1 & cols_non_dd <= len_x_new), ...
    'xtemp: some indices out of bounds in sdp.x.');
assert(length(cols_non_dd) == len_non_ndd, ...
    'xtemp: mismatch between expected and actual length.');

% build selection matrix (1 in each row at the column where the variable lives)
% equivalent to xtemp = jacobian([xlin; xlor; xrot; xpsd], sdp.x)
map_non_dd_x = sparse(1:len_non_ndd, cols_non_dd, 1, len_non_ndd, len_x_new);


% figure out the columns in sdp.g where those non-DD constraints now live
cols_xlin = 1:mlin0;                                            % unchanged
cols_xlor = (mlin0 + m_inserted) + (1:mlor);                    % shifted by m_inserted
cols_xrot = (mlin0 + m_inserted + mlor) + (1:mrot);
cols_xpsd = (mlin0 + m_inserted + mlor + mrot) + (1:mpsd);

cols_non_dd = [cols_xlin, cols_xlor, cols_xrot, cols_xpsd];

% sanity checks for sdp.g
assert(all(cols_non_dd >= 1 & cols_non_dd <= len_g_new), ...
    'map_non_dd_g: some indices out of bounds in sdp.g.');
assert(length(cols_non_dd) == len_non_mdd, ...
    'map_non_dd_g: mismatch between expected and actual length.');
assert(length(unique(cols_non_dd)) == length(cols_non_dd), ...
    'map_non_dd_g: duplicate indices detected in sdp.g mapping.');

% build selection matrix (1 in each row at the column where the variable lives)
% equivalent to gtemp = jacobian([glin; glor; grot; gpsd], sdp.g)
map_non_dd_g = sparse(1:len_non_mdd, cols_non_dd, 1, len_non_mdd, len_g_new);

% empty map
map_dd_eqs_x = sparse(ndd2, len_g_new);
map_dd_eqs_g = sparse(mdd2, len_g_new);

% DD part to deal with the constraints
n_eq_c = dd_index_g.num_eq;       % number of equalities from processing Kc
n_eq_x = dd_index_x.num_eq;       % number of equalities from processing Kx

rows_x = 1:n_eq_x;                % choose first n_eq rows (or any specific ones)
rows_g = 1:n_eq_c;

cols_g = dd_index_g.eq_idx;       % columns from eq_idx
cols_x = dd_index_x.eq_idx;       % columns from eq_idx

% place ones along diagonal
if ndd2~=0
    map_dd_eqs_x(sub2ind([ndd2, len_g_new], rows_x, cols_x)) = 1;
end
if mdd2~=0
    map_dd_eqs_g(sub2ind([mdd2, len_g_new], rows_g, cols_g)) = 1;
end

% construct the combined selection matrix for the dual variables.
% dual variables associated with x
% Non-DD part
lam_x_non_dd = [map_non_dd_x, sparse(len_non_ndd, len_g_new)];

% DD part
comm_Ndd     = blockCommutation(Ndd);
lam_x_dd     = [sparse(ndd2, len_x_new), 0.5 * (speye(ndd2) + comm_Ndd) * map_dd_eqs_x];

% Combine lambda_x components
lam_x = [lam_x_non_dd; lam_x_dd];

% dual variables associated with g
% Non-DD part
lam_g_non_dd = map_non_dd_g;

% DD part
comm_Mdd     = blockCommutation(Mdd);
lam_g_dd     = 0.5 * (speye(mdd2) + comm_Mdd) * map_dd_eqs_g;

% Combine lambda_g components
lam_g = [lam_g_non_dd; lam_g_dd];

% Assemble the full mapping
map.lam = [
    lam_x;
    sparse(len_g_orig, len_x_new), lam_g
];

end
