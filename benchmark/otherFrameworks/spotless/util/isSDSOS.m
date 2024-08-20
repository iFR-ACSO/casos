function [issdsos, Q, phi] = isSDSOS(p)

x = decomp(p);

prog = spotsosprog;
prog = prog.withIndeterminate(x);

[prog,pdsos] = prog.newFreePoly(monomials(x,0:deg(p,x)));
[prog,ind1] = prog.withSDSOS(p);

prog = prog.withPolyEqs(p - pdsos);


% Options
options = spot_sdp_default_options();
options.verbose = 2;
options.solveroptions.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)
options.scale_monomials = false;



sol = prog.minimize(0, @spot_mosek, options);

if (sol.isPrimalFeasible) && (sol.isDualFeasible)
    issdsos = true;
    Q = sol.eval(sol.gramMatrices{ind1});
    phi = sol.gramMonomials{ind1};
else
    issdsos = false;
    Q = [];
    phi = [];
end

end
