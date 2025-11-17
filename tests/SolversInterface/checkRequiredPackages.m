function output = checkRequiredPackages(out_select, packages)
    % Define a list of required packages and their check functions
    % packages = {'mosek'};
    packageChecks = @(pkgName) exist(pkgName, 'file') == 2;

    % Initialize the flag for package availability
    available = true;
    missingPackages = {};

    % Check each package
    for i = 1:numel(packages)
        pkgName = packages{i};
        if ~packageChecks(pkgName)
            available = false;
            missingPackages{end+1} = pkgName;
        end
    end

    % set output
    if out_select==1
        output = available;
    else
        output = missingPackages;
    end
end