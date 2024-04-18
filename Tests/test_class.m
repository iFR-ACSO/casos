clc
clear


%% setup a testsuite for basic math operations
suite_basicMath = testsuite({'test_plus'});

% suite_basicMath = testsuite({'test_plus',...
%                              'test_times',...
%                              'test_power'});

results_suite_basicMath = run(suite_basicMath);


if  all([results_suite_basicMath.Passed])
    disp('Passed all test for basic polynomial math operations.')
else
  
end




