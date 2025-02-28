function [totalTimes,resultTable] = postProcessSolver(varargin)

if isempty(varargin{1})
    error('Provide Solver object for post-processing')
else
    stats = varargin{1}.stats;
end

if nargin >= 2
   plottingFlag = false;
end

% Preallocate arrays for the timing statistics
iterationCount = (1:stats.iterations)'; % Create a column for iteration count
mainSDP = zeros(stats.iterations, 1);
totalBackstepTime = zeros(stats.iterations, 1);
FilterAcceptTime = zeros(stats.iterations, 1);
SuffDecreaseTime = zeros(stats.iterations, 1);
HessApproxTime = zeros(stats.iterations, 1);
SocTime = zeros(stats.iterations, 1);

for k = 1:stats.iterations
    mainSDP(k) = stats.single_iterations{k}.timeStats.mainSDP;
    totalBackstepTime(k) = stats.single_iterations{k}.timeStats.totalBackstepTime;
    FilterAcceptTime(k) = stats.single_iterations{k}.timeStats.FilterAcceptTime;
    SuffDecreaseTime(k) = stats.single_iterations{k}.timeStats.SuffDecreaseTime;
    HessApproxTime(k) = stats.single_iterations{k}.timeStats.HessApproxTime;
    SocTime(k) = stats.single_iterations{k}.timeStats.SocTime;
end

% Create the first table with all timing statistics
resultTable = table(iterationCount, mainSDP, totalBackstepTime, FilterAcceptTime, ...
                    SuffDecreaseTime, HessApproxTime, SocTime, ...
                    'VariableNames', {'iterationCount', 'mainSDP', 'BacksteppingTotal', ...
                                      'FilterAcceptance', 'SufficientDecrease', ...
                                      'HessianApproximation', 'SecondOrderCorrection'});

% Display the first table (optional)
% disp(resultTable);

% Sum the values for each of the timing statistics from the first table
totalMainSDP            = sum(mainSDP);
totalBackstepTime       = sum(totalBackstepTime);
totalFilterAcceptTime   = sum(FilterAcceptTime);
totalSuffDecreaseTime   = sum(SuffDecreaseTime);
totalHessApproxTime     = sum(HessApproxTime);
totalSocTime            = sum(SocTime);

initTime = stats.initTime;
% Create the second table for total times
totalTimes = table(totalMainSDP, totalBackstepTime, totalFilterAcceptTime, ...
                   totalSuffDecreaseTime, totalHessApproxTime, totalSocTime, ...
                   'VariableNames', {'total_mainSDP', 'total_BacksteppingTotal', ...
                                     'total_FilterAcceptance', 'total_SufficientDecrease', ...
                                     'total_HessianApproximation', 'total_SecondOrderCorrection'});


if totalSocTime == 0
    totalSocTime = [];
    explode = [1 1 1 1 1 1 ];
    legendString = {'Initialization','SDP','Filter Acceptance','Sufficient Decrease', 'Hessian Approxmation','Other'};
else
    explode = [1 1 1 1 1 1 1];
    legendString = {'Initialization','SDP','Filter Acceptance','Sufficient Decrease', 'Hessian Approxmation', 'Second-order-Correction','Other'};
end
timesSub = sum([initTime,totalMainSDP,totalFilterAcceptTime,totalSuffDecreaseTime,totalHessApproxTime,totalSocTime]);

% verbosity, simple function evaluations, function calls, ...
other = stats.totalSolveTime -  timesSub;

normalizedTimes = [initTime,totalMainSDP,totalFilterAcceptTime,totalSuffDecreaseTime,totalHessApproxTime,totalSocTime,other]./sum([initTime,totalMainSDP,totalFilterAcceptTime,totalSuffDecreaseTime,totalHessApproxTime,totalSocTime,other]);
figure('Name','Timing Statistics')
pie(normalizedTimes,explode)
legend(legendString,'Location','westoutside')

end
