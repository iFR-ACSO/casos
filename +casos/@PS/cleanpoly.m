function B = cleanpoly(P,tol,deg)

%% check inputs
if nargin<3
    deg = [];
end


assert(~is_symbolic(P.coeffs), 'Coefficients must be numerical values!')
assert(isempty(tol) ||((isa(tol,'double') && isscalar(tol) && tol>=0)), 'Tolerance must be a positive scalar or is empty!')

p = casos.PS(P);

%% remove numeric coefficients that are below tol
if ~isempty(tol)
  
        coeffs = p.coeffs;

        % polynomial coefficients are casadi SX; we need numerical values
        idx = find(abs(full(casadi.DM(coeffs))) < tol);

        if ~isempty(idx)
            
            % remove coefficients that are below tolerance
            coeff      = full(casadi.DM(coeffs));
            coeff(idx) = 0;

            % new polynomial
            B        = casos.PS;

            % remove zero terms
            [coeffs,degmat,indets] = removeZero(sparsify(casadi.SX(coeff)),p.degmat,p.indets);

            B.coeffs = casadi.SX(coeffs); % polynomial coefficients need to be casadi SX
            B.degmat = degmat; 
            B.indets = indets;
            B.matdim = p.matdim;

        else
            % in case all terms are above tolerance
            B = p;

        end
end

%% find terms we want to retain 

if ~isempty(deg)
    
    % find all parts that shall be retained with certain degrees
    if isa(deg,'double') && size(deg,2) == 2 && ...
        all(floor(deg)==ceil(deg)) && all(deg>=0)

        
        pdeg = sum(p.degmat,2);

        % find all terms we want to retain
        idx = arrayfun(@(x) find(pdeg == x), unique(deg), 'UniformOutput', false);
        idx = vertcat(idx{:});
        
        % find all terms we don't want to retain
        idx = setdiff(1:length(pdeg),idx);
        
        % remove the coefficients
        coeffs        =  p.coeffs;
        coeffs(idx,:) = 0;

        % new polynomial
        B        = casos.PS;

       % remove zero terms
       [coeffs,degmat,indets] = removeZero(sparsify(casadi.SX(coeffs)),p.degmat,p.indets);

       B.coeffs = coeffs; % polynomial coefficients need to be casadi SX
       B.degmat = degmat; 
       B.indets = indets;
       B.matdim = p.matdim;
    
    % find all parts that shall be retained with certain degrees of a
    % variable
    elseif iscell(deg)  && size(deg,2)== 2

        % deg is a cell array of pvars and vectors of non-neg integers
        coeffs = p.coeffs;

        for i1=1:size(deg,1)
            
            % check if variable is part of indeterminates
            var = deg{i1,1};

            vidx = find(strcmp(p.indets,var));
            
            if ~isempty(vidx)
                
                degi = deg{i1,2};
                
                % check which degrees shall be retained
                if ~(isa(degi,'double') && size(degi,2) == 2 && ...
                      all(floor(degi) == ceil(degi)) && all(degi>=0) )

                    error(['Second column of cell array deg must '...
                           'contain a vector of non-negative integers.']);
                end
                
                % get entries in degmat of the desired indeterminates we
                % want to retain
                pdeg = p.degmat(:,vidx);


                idx = arrayfun(@(x) find(pdeg == x), unique(degi), 'UniformOutput', false);
                idx = vertcat(idx{:});

                % get entries we don't want to retain
                idx = setdiff(1:length(pdeg),idx);
                coeffs(idx,:) = 0;

            end
        end

        % new polynomial
        B        = casos.PS;

       % remove zero terms
       [coeffs,degmat,indets] = removeZero(sparsify(casadi.SX(coeffs)),p.degmat,p.indets);

       B.coeffs = coeffs; % polynomial coefficients need to be casadi SX
       B.degmat = degmat; 
       B.indets = indets;
       B.matdim = p.matdim;

       


    else

        error('deg must be a vector of non-negative integers or an Nx2 cell array');

    end

    
    
end


