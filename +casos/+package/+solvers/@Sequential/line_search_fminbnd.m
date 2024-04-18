function dopt= line_search_fminbnd(obj, xk, xk1, sol)
 
   
    dopt = fminbnd(@(d) double(obj.Merit(xk.*(1-d) + xk1.*d, sol{5})), 0.1, 1);
  
end