function [sdp,args,map,opts] = sdd_reduce(obj, sdp, opts, args)
% Reduce a SDD cone program to SOCP 

% check cones
check_cone(obj.get_cones,opts.Kx,'lin');
check_cone(obj.get_cones,opts.Kc,'lin');
check_cone(obj.get_cones,opts.Kx,'sdd');
check_cone(obj.get_cones,opts.Kc,'sdd');

% get dimensions of cones in the program decision variables
Nlin = get_dimension(obj.get_cones,opts.Kx,'lin');
Nsdd = get_dimension(obj.get_cones,opts.Kx,'sdd');
% (unused, but collected for completeness)
Nlor = get_dimension(obj.get_cones,opts.Kx,'lor');
Nrot = get_dimension(obj.get_cones,opts.Kx,'rot');
Npsd = get_dimension(obj.get_cones,opts.Kx,'psd');
Ndd  = get_dimension(obj.get_cones,opts.Kx,'dd');


% get dimensions of cones in the program constraints
Mlin = get_dimension(obj.get_cones,opts.Kc,'lin');
Msdd = get_dimension(obj.get_cones,opts.Kc,'sdd');
% (unused, but collected for completeness)
Mlor = get_dimension(obj.get_cones,opts.Kc,'lor');
Mrot = get_dimension(obj.get_cones,opts.Kc,'rot');
Mpsd = get_dimension(obj.get_cones,opts.Kc,'psd');
Mdd  = get_dimension(obj.get_cones,opts.Kc,'dd');

% temporarily save sdp.x (original)
x_original = sdp.x;
g_original = sdp.g;

[xlin, xlor, xrot, xpsd, xdd, xsdd] = separate(x_original, [sum(Nlin) sum(Nlor) sum(Nrot) sum(Npsd) sum(Ndd.^2) sum(Nsdd.^2)]);

% initialize SDD index structs in case the fields don't exist
sdd_index_g.num_eq = 0;
sdd_index_x.num_eq = 0;

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
sdp.x = [sdp.x(1:Nlin); added_vars; sdp.x(Nlin+1:end)];

% add the variables M to sdp.x
Nlin = Nlin + sum(Nsdd.^2)+length(added_vars);
Mlin = Mlin + sdd_index_x.num_eq+sdd_index_g.num_eq; 
Mpsd = [Mpsd, repmat(2, 1, length(added_vars)/3)];


% update lbx and ubx 
args.dd_lbx = [args.dd_lbx; -inf(sum(Nsdd.^2) + length(added_vars),1)];
args.dd_ubx = [args.dd_ubx;  inf(sum(Nsdd.^2) + length(added_vars),1)];

% update linear variables and constraints
opts.Kx = setfield(opts.Kx,'lin',Nlin);
opts.Kc = setfield(opts.Kc,'lin',Mlin);

% update linear variables and constraints
opts.Kx = setfield(opts.Kx,'lin', Nlin);
opts.Kc = setfield(opts.Kc,'lin', Mlin);
opts.Kc = setfield(opts.Kc,'psd', Mpsd);

% map from new sdp.x to old (only the non-sdd variables)
map.x = jacobian(x_original, sdp.x);

% map from new sdp.g to old
map.g = [speye(length(g_original)), sparse(length(g_original), length(sdp.g)-length(g_original))];

% map from new lam_* to original lam
xtemp = jacobian([xlin; xlor; xrot; xpsd; xdd], sdp.x);
G = sparse(length(xsdd), size(map.g,2));
n_eq = sdd_index_x.num_eq;              % number of equalities
rows = 1:n_eq;                          % choose first n_eq rows (or any specific ones)
cols = sdd_index_x.eq_idx;              % columns from eq_idx

% place ones along diagonal
G(sub2ind(size(G), rows, cols)) = 1;

map.lam_x = [xtemp, sparse(size(xtemp, 1), size(map.g,2));
             zeros(length(xsdd), size(xtemp,2)), G];

map.lam_a = [sparse(size(map.g,1), size(xtemp,2)), map.g];

end


function varargout = separate(A,varargin)
% Separate array into subarrays.
    varargout = mat2cell(A,varargin{:});
end




