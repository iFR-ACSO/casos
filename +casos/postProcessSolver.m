function [totalTimes,resultTable] = postProcessSolver(varargin)

if isempty(varargin{1})
    error('Provide Solver object for post-processing')
else
    stats = varargin{1}.stats;
end

if nargin > 2
   plottingFlag = false;
else
    plottingFlag = varargin{2};
end

% Preallocate arrays for the timing statistics
iterationCount = (1:stats.iterations)'; % Create a column for iteration count
mainSDP = zeros(stats.iterations, 1);
totalBackstepTime = zeros(stats.iterations, 1);
FilterAcceptTime = zeros(stats.iterations, 1);
SuffDecreaseTime = zeros(stats.iterations, 1);
HessApproxTime = zeros(stats.iterations, 1);
SocTime = zeros(stats.iterations, 1);
FeasResTime = zeros(stats.iterations, 1);
for k = 1:stats.iterations
    mainSDP(k) = stats.single_iterations{k}.timeStats.mainSDP;
    totalBackstepTime(k) = stats.single_iterations{k}.timeStats.totalBackstepTime;
    FilterAcceptTime(k) = stats.single_iterations{k}.timeStats.FilterAcceptTime;
    SuffDecreaseTime(k) = stats.single_iterations{k}.timeStats.SuffDecreaseTime;
    HessApproxTime(k) = stats.single_iterations{k}.timeStats.HessApproxTime;
    SocTime(k) = stats.single_iterations{k}.timeStats.SocTime;
    FeasResTime(k) = stats.single_iterations{k}.timeStats.feasResTime;
end

% Create the first table with all timing statistics
resultTable = table(iterationCount, mainSDP, totalBackstepTime, FilterAcceptTime, ...
                    SuffDecreaseTime, HessApproxTime, SocTime,FeasResTime, ...
                    'VariableNames', {'iterationCount', 'main-SDP', 'Backstepping-Total', ...
                                      'Filter-Acceptance', 'Sufficient-Decrease', ...
                                      'Hessian-Approximation', 'Second-Order-Correction', 'Feasibility-Restoration'});


% Sum the values for each of the timing statistics from the first table
totalMainSDP            = sum(mainSDP);
totalBackstepTime       = sum(totalBackstepTime);
totalFilterAcceptTime   = sum(FilterAcceptTime);
totalSuffDecreaseTime   = sum(SuffDecreaseTime);
totalHessApproxTime     = sum(HessApproxTime);
totalSocTime            = sum(SocTime);
FeasResTime             = sum(FeasResTime);

% Create the second table for total times
totalTimes = table(totalMainSDP, totalBackstepTime, totalFilterAcceptTime, ...
                   totalSuffDecreaseTime, totalHessApproxTime, totalSocTime,FeasResTime , ...
                   'VariableNames', {'total_main-SDP', 'total_Backstepping', ...
                                     'total_Filter-Acceptance', 'total_Sufficient-Decrease', ...
                                     'total_Hessian-Approximation', 'total_Second-Order-Correction','total_Feasibility-Restoration'});


if plottingFlag
   initTime = stats.initTime;

% second-order correction and/or feas. restoration are not necessarily
% called
if totalSocTime == 0
    totalSocTime = [];
    explode = [1 1 1 1 1 1 1];
    legendString = {'Initialization','SDP','Filter Acceptance','Sufficient Decrease', 'Hessian Approxmation','Feasbility-Restoration','Other'};
elseif FeasResTime == 0
    FeasResTime = [];
    explode = [1 1 1 1 1 1 1];
    legendString = {'Initialization','SDP','Filter Acceptance','Sufficient Decrease', 'Hessian Approxmation', 'Second-order-Correction','Other'};
elseif  FeasResTime == 0 && totalSocTime == 0
    FeasResTime = [];
    totalSocTime = [];

    explode = [1 1 1 1 1 1];
    legendString = {'Initialization','SDP','Filter Acceptance','Sufficient Decrease', 'Hessian Approxmation','Other'};
else
    explode = [1 1 1 1 1 1 1 1];
    legendString = {'Initialization','SDP','Filter Acceptance','Sufficient Decrease', 'Hessian Approxmation', 'Second-order-Correction','Feasbility-Restoration','Other'};
end

timesSub = [initTime,totalMainSDP,totalFilterAcceptTime,totalSuffDecreaseTime,totalHessApproxTime,totalSocTime, FeasResTime];

% verbosity, simple function evaluations,termination criteria, function calls, ...
other = stats.totalSolveTime -  sum(timesSub);

normalizedTimes = [timesSub,other]./sum([timesSub,other]);

figure('Name','Timing Statistics')
pie(normalizedTimes,explode)
legend(legendString,'Location','westoutside')
end

end
