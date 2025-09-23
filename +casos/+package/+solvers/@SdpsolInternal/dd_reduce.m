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
% (unused, but collected for completeness)
Nsdd = get_dimension(obj.get_cones,opts.Kx,'sdd'); nsdd2  = sum(Nsdd.^2);

% get dimensions of cones in the program constraints
Mlin = get_dimension(obj.get_cones,opts.Kc,'lin'); mlin0 = sum(Mlin);
Mdd  = get_dimension(obj.get_cones,opts.Kc,'dd');  mdd2  = sum(Mdd.^2);
% (unused, but collected for completeness)
Mlor = get_dimension(obj.get_cones,opts.Kc,'lor'); mlor  = sum(Mlor);
Mrot = get_dimension(obj.get_cones,opts.Kc,'rot'); mrot  = sum(Mrot);
Mpsd = get_dimension(obj.get_cones,opts.Kc,'psd'); mpsd  = sum(Mpsd.^2);
Msdd = get_dimension(obj.get_cones,opts.Kc,'sdd'); msdd2 = sum(Msdd.^2);

% save sizes of sdp.x and sdp.g (original)
len_x_orig = length(sdp.x);     % original decision variable length
len_g_orig = length(sdp.g);     % original constraint length

% create zero initial 
num_nlin_x = 0;     % in decision variables
num_nlin_g = 0;     % in constraints

% verify DD cones in the constraints and create slack DD variables
M_g = cell(length(Mdd),1);
if isfield(opts.Kc,'dd')
    [sdp,args,M_g,num_nlin_g,dd_index_g,opts] = obj.replaceDDcones(sdp, Mdd, Mlin, args, opts, 'g');
else
    dd_index_g.num_eq = 0;
    dd_index_g.eq_idx = [];
end

% verify DD cones in decision variables
M_x = cell(length(Ndd),1);
if isfield(opts.Kx, 'dd')
    [sdp,args,M_x,num_nlin_x,dd_index_x, opts] = obj.replaceDDcones(sdp, Ndd, Mlin, args, opts, 'x');
else
    dd_index_x.num_eq = 0;
    dd_index_x.eq_idx = [];
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
cols_xlor = (mlin0 + num_nlin_x + num_nlin_g) + (1:mlor);                    % shifted by n_inserted
cols_xrot = (mlin0 + num_nlin_x + num_nlin_g + mlor) + (1:mrot);
cols_xpsd = (mlin0 + num_nlin_x + num_nlin_g + mlor + mrot) + (1:mpsd);

cols_non_dd = [cols_xlin, cols_xlor, cols_xrot, cols_xpsd];

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
% treat lam_x
map.lam = [ ...
    map_non_dd_x,                sparse(len_non_ndd, len_g_new);     % lam_x (non-DD part)
    sparse(ndd2, len_x_new),     map_dd_eqs_x ...                   % lam_x (DD part)
];

map.lam = [map.lam;
           sparse(len_g_orig, len_x_new), [ map_non_dd_g;                                               % lam_g (original constraints)
                                            0.5*(speye(mdd2)+blockCommutation(Mdd))*map_dd_eqs_g]];     % map to the duals of DD 

end


function K = blockCommutation(Mdd)
% Mdd = [n1, n2, ..., nk] where each ni is a block size
% Build commutation matrices for each block
Kblocks = cellfun(@commutationMatrix, num2cell(Mdd), 'UniformOutput', false);

% Assemble block diagonal matrix
K = blkdiag(Kblocks{:});

% Ensure sparse format
K = sparse(K);
end

function K = commutationMatrix(n)
% commutationMatrix: Construct commutation matrix K_{n,n}.
%   K = commutationMatrix(n) returns the sparse commutation matrix of size
%   n^2-by-n^2 such that:
%       K * vec(X) = vec(X.')

% Preallocate indices
rows = zeros(n^2,1);
cols = zeros(n^2,1);

% Compute permutation mapping
idx = 1;
for i = 1:n
    for j = 1:n
        % vec(X) index for (i,j)
        col = (j-1)*n + i;
        % vec(X.') index for (i,j)
        row = (i-1)*n + j;

        rows(idx) = row;
        cols(idx) = col;
        idx = idx + 1;
    end
end

% Build sparse permutation matrix
K = sparse(rows, cols, 1, n^2, n^2);
end


