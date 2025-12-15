function [sdp,args,map,opts] = sdd_reduce(obj,sdp,opts,args)
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
Npsd = get_dimension(obj.get_cones,opts.Kx,'psd'); npsd   = sum(Npsd.^2);
% (unused, but collected for completeness)
Ndd  = get_dimension(obj.get_cones,opts.Kx,'dd');  ndd2   = sum(Ndd.^2);

% get dimensions of cones in the program constraints
Mlin = get_dimension(obj.get_cones,opts.Kc,'lin'); mlin0 = sum(Mlin);
Msdd = get_dimension(obj.get_cones,opts.Kc,'sdd'); msdd2 = sum(Msdd.^2);
Mpsd = get_dimension(obj.get_cones,opts.Kc,'psd'); mpsd  = sum(Mpsd.^2);
Mdd  = get_dimension(obj.get_cones,opts.Kc,'dd');  mdd2  = sum(Mdd.^2);
% (unused, but collected for completeness)
Mlor = get_dimension(obj.get_cones,opts.Kc,'lor'); mlor  = sum(Mlor);
Mrot = get_dimension(obj.get_cones,opts.Kc,'rot'); mrot  = sum(Mrot);

% save sizes of sdp.x and sdp.g (original)
len_x_orig = size(sdp.x, 1);     % original decision variable length
len_g_orig = size(sdp.g, 1);     % original constraint length

% initialize SDD index structs in case the fields don't exist
sdd_index_g.num_eq = 0;         % in constraints
sdd_index_x.num_eq = 0;         % in decision variables

% verify SDD cones in the constraints and create slack SDD variables
M_g = cell(numel(Msdd),1);
if msdd2 > 0
    [sdp,args,M_g,~,sdd_index_g,opts] = replaceSDDcones(sdp, Msdd, Mlin, args, opts, 'g');
else
    sdd_index_g.num_eq = 0;
    sdd_index_g.eq_idx = [];
end

% verify SDD cones in decision variables
M_x = cell(length(Nsdd),1);
if nsdd2 > 0
    [sdp,args,M_x,~,sdd_index_x, opts] = replaceSDDcones(sdp, Nsdd, Mlin, args, opts, 'x');
else
    sdd_index_x.num_eq = 0;
    sdd_index_x.eq_idx = [];
end

% update decision variables
added_vars = [vertcat(M_g{:}); vertcat(M_x{:})];
n_inserted = numel(added_vars);

sddvar  = sdp.x(nlin0+nlor+nrot+npsd+ndd2+1:end);
restvar = sdp.x(nlin0+nlor+nrot+1:nlin0+nlor+nrot+npsd+ndd2);
sdp.x = [sdp.x(1:Nlin); added_vars; sddvar; restvar];

n_eq_c = sdd_index_g.num_eq;   % number of equalities from processing Kc

% add the variables M to sdp.x
Nlin = Nlin + nsdd2+n_inserted;
Mlin = Mlin + sdd_index_x.num_eq+sdd_index_g.num_eq; 

% combine both lists (processed separately but same logic)
allSdd = {Msdd, Nsdd};

% compute total length needed
totalLen = 0;
for k = 1:numel(allSdd)
    v = allSdd{k};
    totalLen = totalLen + sum(max(1, v .* (v - 1) / 2));
end

% preallocate
Mpsd = zeros(1, totalLen);
pos  = 1;

for k = 1:numel(allSdd)
    v = allSdd{k}(end:-1:1);      % reverse once

    numPairs = v .* (v - 1) / 2;
    counts   = max(numPairs, 1); % how many entries each cone contributes

    % total length needed for this cell
    L = sum(counts);

    % prefill with 2's
    Mpsd(pos:pos+L-1) = 2;

    % positions where numPairs == 0 should be 1 instead of 2
    zeroIdx = find(numPairs == 0);
    if ~isempty(zeroIdx)
        offsets = cumsum([0; counts(1:end-1)]);
        Mpsd(pos + offsets(zeroIdx)) = 1;
    end

    pos = pos + L;
end

% update lbx and ubx 
args.dd_lbx = [args.dd_lbx; -inf(nsdd2 + n_inserted,1)];
args.dd_ubx = [args.dd_ubx;  inf(nsdd2 + n_inserted,1)];

% update cones
opts.Kx.lin = Nlin;
opts.Kc.lin = Mlin;
opts.Kc.psd = Mpsd;

% get sizes of NEW sdp.x and sdp.g
len_x_new  = numel(sdp.x);
len_g_new  = numel(sdp.g);

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
len_non_sdd_x   = nlin0 + nlor + nrot + npsd;
len_non_sdd_g = mlin0 + mlor + mrot + mpsd + mdd2;

% figure out the columns in sdp.x where non-SDD variables now live
cols_xlin = 1:nlin0;                                          % unchanged
cols_xlor = (nlin0 + n_inserted) + (1:nlor);                  % shifted by inserted vars
cols_xrot = (nlin0 + n_inserted + nlor) + (1:nrot);
cols_xpsd = (nlin0 + n_inserted + nlor + nrot) + (1:npsd);

cols_non_sdd = [cols_xlin, cols_xlor, cols_xrot, cols_xpsd];

% sanity checks
assert(all(cols_non_sdd >= 1 & cols_non_sdd <= len_x_new), ...
    'xtemp: some indices out of bounds in sdp.x.');
assert(length(cols_non_sdd) == len_non_sdd_x, ...
    'xtemp: mismatch between expected and actual length.');

% selection matrix for non-SDD variables
% equivalent to xtemp = jacobian([xlin; xlor; xrot; xpsd; xdd], sdp.x);
map_non_sdd_x = sparse(1:len_non_sdd_x, cols_non_sdd, 1, len_non_sdd_x, len_x_new);

% figure out indices where these now live in sdp.g
cols_glin = 1:mlin0;                        % linear constraints unchanged
cols_glor = (mlin0 + n_eq_c) + (1:mlor);    % shifted by inserted equalities
cols_grot = (mlin0 + n_eq_c + mlor) + (1:mrot);
cols_gpsd = (mlin0 + n_eq_c + mlor + mrot) + (1:mpsd);
cols_gdd  = (mlin0 + n_eq_c + mlor + mrot + mpsd) + (1:mdd2);

cols_non_sdd_g = [cols_glin, cols_glor, cols_grot, cols_gpsd, cols_gdd];

% sanity check
assert(length(cols_non_sdd_g) == len_non_sdd_g, ...
    'Mismatch between expected and actual non-SDD g length.');

% selection matrix for non-SDD constraints
map_non_sdd_g = sparse(1:len_non_sdd_g, cols_non_sdd_g, 1, len_non_sdd_g, len_g_new);


% empty map
map_sdd_eqs_x = sparse(nsdd2, len_g_new);
map_sdd_eqs_g = sparse(msdd2, len_g_new);

% SDD part to deal with constraints
n_eq_c = sdd_index_g.num_eq;       % number of equalities from processing Kc
n_eq_x = sdd_index_x.num_eq;       % number of equalities from processing Kx

rows_x = 1:n_eq_x;                % choose first n_eq rows (or any specific ones)
rows_g = 1:n_eq_c;

cols_g = sdd_index_g.eq_idx;       % columns from eq_idx
cols_x = sdd_index_x.eq_idx;       % columns from eq_idx

% place ones along diagonal
if nsdd2~=0
    map_sdd_eqs_x(sub2ind([nsdd2, len_g_new], rows_x, cols_x)) = 1;
end
if msdd2~=0
    map_sdd_eqs_g(sub2ind([msdd2, len_g_new], rows_g, cols_g)) = 1;
end

% construct the combined selection matrix for the dual variables.
% dual variables associated with x
% Non-SDD part
lam_x_non_sdd = [map_non_sdd_x, sparse(len_non_sdd_x, len_g_new)];

% SDD part
comm_Nsdd     = blockCommutation(Nsdd);
lam_x_sdd     = [sparse(nsdd2, len_x_new), 0.5 * (speye(nsdd2) + comm_Nsdd) * map_sdd_eqs_x];

% Combine lambda_x components
lam_x = [lam_x_non_sdd; lam_x_sdd];

% dual variables associated with g
% Non-SDD part
lam_g_non_sdd = map_non_sdd_g;

% SDD part
comm_Msdd     = blockCommutation(Msdd);
lam_g_sdd     = 0.5 * (speye(msdd2) + comm_Msdd) * map_sdd_eqs_g;

% Combine lambda_g components
lam_g = [lam_g_non_sdd; lam_g_sdd];

% Assemble the full mapping
map.lam = [
    lam_x;
    sparse(len_g_orig, len_x_new), lam_g
];

end
