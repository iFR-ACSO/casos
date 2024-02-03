function b = int(a,x,L,U)

if nargin==1
    x = a.indeterminates.str;
    L = []; 
    U = [];
elseif nargin==2
    if isa(x,'double')
        % B = int(A,[L,U])
        L = x(1);
        U = x(2);
        x = a.indeterminates.str;
    else
        % B = int(A,X)
        L = []; U = [];
    end
elseif nargin==3
    if isa(x,'double')
        % B = int(A,L,U)
        U = L;
        L = x;
        x = a.indeterminates.str;
    else
        % B = int(A,X,[L U])
        U = L(2);
        L = L(1);
    end
elseif nargin~=4
    error(['Invalid syntax for the "int" command. ' ...
        'Type "help int" for more information.'])
end

if length(x)==1
    x = x.indeterminates;
elseif ~ischar(x.indeterminates.str)
    error('X must be a single polynomial variable or a string');
end

a = casos.PS(a);
b = casos.PS();


% Get polynomial info about A
acoef   = a.coeffs;
adegmat = a.degmat;
avar    = a.indeterminates;
adim    = a.matdim;
Naterms = size(acoef,1);

% Find variable we are differentiating with respect to.
varnumb = find( strcmp(a.indeterminates.str, x.indeterminates.str) );

% Perform indefinite integral
bdim = adim;

if isempty(varnumb)
    % int a(y) dx = a(y)*x
    bcoef = acoef;
    bdegmat = [adegmat ones(Naterms,1)];
    bvar = [avar; x];
else
    % int a(y)*x^n dx = a(y)*x^(n+1) / (n+1)
    bdegmat = adegmat;
    nplus1 = bdegmat(:,varnumb)+1;
    bdegmat(:,varnumb) = nplus1;
    
    bcoef = acoef;
    bcoef = lrscale( bcoef, 1./nplus1 , []);
    
    bvar = avar;
end
    
b.coeffs = bcoef;
b.degmat = bdegmat;
b.indets = str(bvar);
b.matdim = bdim ;

% Evaluate definite integral
if ~isempty(L)
    bU = subs(b,x,U);
    bL = subs(b,x,L);
    b = bU-bL;
end

end % end of function