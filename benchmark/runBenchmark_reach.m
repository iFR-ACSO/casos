% script that runs the GTM 4D ROA problem with different toolboxes

clc
clear
close all

% add paths to subfolder
addpath('./casos_ex/')
addpath('./yalmip_ex/')
addpath('./sostools_ex/')
addpath('./sosopt_ex/')
addpath('./spotless_ex/')


[gval_sopt_GTM,solverTime_total_sopt_GTM,buildTime_sopt_GTM] = reachEstGTM_benchSosopt();

[gval_c_GTM,solverTime_total_c_GTM,buildTime_c_GTM]= reachEstGTM_benchCasos();

