function  pcontour(p, g, opts)
% Plots contour of p(x) with value g.
% The contours are generated numerically by evaluating p on a grid of 
% values in each direction. 
% Option 'auto' in opts allows the method to obtain a fitted bound on the
% level-set

var = get_indets(p);    % obtain indeterminates 
dim = var.length;       % dimension 

% TODO: add an option to slice in order to become either 2 or 3 dimensional
assert(dim<3,'p must be a polynomial in 2 or 3 variables')

if nargin == 1
    g = 1;
end

if nargin < 3
    opts = struct( ...
        'domain',   repmat([-1 1], 1, dim) ,...
        'linespec', 'k' ,...
        'npts',     repmat(50, 1, dim),...
        'auto',     1 ...
        );
end

% automatically obtain domain
if opts.auto == 1
    % auxialiary variables
    b = casos.PS.sym('b');
    sign = casos.PS.sym('sign');
    m = casos.PS.sym('m',[dim,1]);
    s1 = casos.PS.sym('s1', monomials(var,0:2), 'gram');
    
    % define SOS feasibility
    sos = struct('x',[b; s1],'g',+s1*(p-g)+ sign*(sign*b-m'*var), 'f', b, 'p', [sign; m]);

    % states + constraint are SOS cones
    opts_opt.Kx.l = 1;
    opts_opt.Kx.s = 1;
    opts_opt.Kc.s = 1;
    % ignore infeasibility
    opts_opt.error_on_fail = false;

    % solve by relaxation to SDP
    S = casos.sossol('S','sedumi',sos,opts_opt);
    
    dir = kron(eye(dim), [1 1]);
    for i=1:length(opts.domain)
        sol = S('p', [opts.domain(i); dir(:,i)]);
        opts.domain(i) = opts.domain(i)*sol.x(1);
    end

end
    
% discretize the domain
dg = cell(dim,1);
for i=1:dim
    dg{i} = linspace(opts.domain(2*i-1),opts.domain(2*i),opts.npts(i));
end

% pre-allocate memory for pgrid
pgrid = zeros(opts.npts);
pgrid = pgrid(:);

% convert polynomial to DM function
p = to_function(p);

if dim == 2 %2D
    [dg{1},dg{2}] = meshgrid(dg{1},dg{2});
    xg_c = dg{1}(:);
    yg_c = dg{2}(:);
    for i=1:length(pgrid)
        pgrid(i) = full(p([xg_c(i); yg_c(i)]));
    end
    pgrid = reshape(pgrid, opts.npts);
    contour(dg{1},dg{2},pgrid,[g g],opts.linespec);  
    axis(opts.domain)
    xlabel('x')
    ylabel('y')
else    % 3D
    [dg{1},dg{2},dg{3}] = meshgrid(dg{1},dg{2}, dg{3});
    xg_c = dg{1}(:);
    yg_c = dg{2}(:);
    zg_c = dg{3}(:);
    for i=1:length(pgrid)
        pgrid(i) = full(p([xg_c(i); yg_c(i); zg_c(i)]));
    end
    pgrid = reshape(pgrid, opts.npts);
    isosurface(xg,yg,zg,pgrid,1);    
    axis(opts.domain)
    xlabel('x')
    ylabel('y')
    zlabel('z')
end

end


