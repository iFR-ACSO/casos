function dopt= line_search_fminbnd(obj, p0, xk, xk1, dual_plus )
 
   
    dopt = fminbnd(@(d) double(obj.Merit(xk.*(1-d) + xk1.*d, dual_plus,p0 )), 0.1, 1);
  
end