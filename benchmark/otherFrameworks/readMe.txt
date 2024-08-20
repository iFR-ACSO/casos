Changes in SOSOPT/GSOPT

- Within bisection we store ALL iteration info structs.
- Add tic-toc time measure in sedumi2mosek.m and extended info struct


SOSTOOLS
- Add tic-toc time measure in sedumi2mosek.m and extended info struct
- commented out several parts with disp(); could not be turned off by parameter

SPOTless
-augmented Mosek STATUS switch case with "otherwise"; in some cases an error occured because the case was not defined (probably newer Mosek versions have more outputs)

