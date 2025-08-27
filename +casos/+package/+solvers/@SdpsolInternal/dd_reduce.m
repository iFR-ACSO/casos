function [sdp,args,map,opts] = dd_reduce(obj, sdp, opts, args)
% Reduce a DD cone program to LP 

% check cones
check_cone(obj.get_cones,opts.Kx,'lin');
check_cone(obj.get_cones,opts.Kc,'lin');
check_cone(obj.get_cones,opts.Kx,'dd');
check_cone(obj.get_cones,opts.Kc,'dd');

% get dimensions of cones in the program decision variables
Nlin = get_dimension(obj.get_cones,opts.Kx,'lin');
Ndd  = get_dimension(obj.get_cones,opts.Kx,'dd');
% (unused, but collected for completeness)
Nlor = get_dimension(obj.get_cones,opts.Kx,'lor');
Nrot = get_dimension(obj.get_cones,opts.Kx,'rot');
Npsd = get_dimension(obj.get_cones,opts.Kx,'psd');
Nsdd = get_dimension(obj.get_cones,opts.Kx,'sdd');

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
num_eq_g = 0;
num_eq_x = 0;
num_ineq_g = 0;
num_ineq_x = 0;

% verify DD cones in the constraints and create slack DD variables
M_g = cell(length(Mdd),1);
if isfield(opts.Kc,'dd')
    [sdp,args,M_g,num_eq_x,num_ineq_g,opts] = obj.replaceDDcones(sdp, Mdd, Mlin, args, opts, 'g');
end

% verify DD cones in decision variables
M_x = cell(length(Ndd),1);
if isfield(opts.Kx, 'dd')
    [sdp,args,M_x,num_eq_x,num_ineq_x, opts] = obj.replaceDDcones(sdp, Ndd, Mlin, args, opts, 'x');
end

% update decision variables
added_vars = [vertcat(M_g{:}); vertcat(M_x{:})];
sdp.x = [sdp.x(1:Nlin); added_vars; sdp.x(Nlin+1:end)];

% add linear cones to constraints
Mlin = Mlin + num_ineq_x + num_eq_x + num_eq_g + num_ineq_g;

% add the variables M to sdp.x
n_added = sum(Ndd.^2) + length(added_vars);
Nlin = Nlin + n_added;

args.dd_lbx = [args.dd_lbx; -inf(n_added,1)];
args.dd_ubx = [args.dd_ubx;  inf(n_added,1)];

% update linear variables and constraints
opts.Kx = setfield(opts.Kx,'lin',Nlin);
opts.Kc = setfield(opts.Kc,'lin',Mlin);

% map from new sdp.x to old
map.x = jacobian(x_original, sdp.x);

% map from new sdp.g to old
map.g = [speye(length(g_original)), sparse(length(g_original), length(sdp.g)-length(g_original))];

% map from new lam_* to original lam
map.lam_x = [map.x, sparse( size(map.x, 1), size(map.g,2))];
map.lam_a = [sparse(size(map.g,1), size(map.x,2)), map.g];

end






