function [Z,K,z,Mp,Md] = grambasis(S,I,newton_solver)
% Return Gram basis of polynomial vector.

if nargin < 3
    % no Newton simplification
    newton_solver = '';
end

if nargin < 2 || isempty(I)
    lp = numel(S);
    I = true(lp,1);
    idx = 1:lp;
else
    lp = nnz(I);
    idx = find(I);
end

% get logical maps for degrees and indeterminate variables
% -> Ldeg(i,j) is true iff p(i) has terms of degree(j)
% -> Lvar(i,j) is true iff p(i) has terms in indets(j)
[degree,Ldeg] = get_degree(S,I);
[indets,Lvar] = get_indets(S,I);
% remove non-indexed logicals
Ldeg(~I,:) = []; Lvar(~I,:) = [];
% detect degrees and variables
Id = any(Ldeg,1); Iv = any(Lvar,1);
% remove unused degrees or variables
Ldeg(:,~Id) = []; Lvar(:,~Iv) = [];

% min and max degree of Gram basis vector
mndg = floor(min(degree)/2);
mxdg =  ceil(max(degree)/2);

% ensure min and max degree are even
[degree,ii] = unique([degree,(2*mndg):(2*mxdg)]);
ld = length(degree); tmp = [Ldeg false(lp,ld)];
Ldeg = tmp(:,ii);

% split into even and odd monomials
Ldeg_e = Ldeg(:,mod(degree,2)==0);
Ldeg_o = Ldeg(:,mod(degree,2) >0);
% add even degrees if necessary
Ldeg_e(:,1:end-1) = Ldeg_e(:,1:end-1) | Ldeg_o;
Ldeg_e(:,2:end) = Ldeg_e(:,2:end) | Ldeg_o;

% get half-degree vector of monomials
% TODO: compute degree matrix directly (avoid unnecessary checks)
z = to_vector(casos.Sparsity.scalar(indets,mndg:mxdg));
% compute logical map for half-degree monomials
% -> Lz_deg(i,j) is true iff z(i) is of degree(j)/2
% -> Lz_var(i,j) is true iff z(i) includes indets(j)
[~,Lz_deg] = get_degree(z);
[~,Lz_var] = get_indets(z);

% build logical map for half-degree monomials
% -> Lz(i,j) is true iff Gram form of p(i) includes z(j)
% that is, z(j) is of a degree in the Gram basis for p(i)
% AND z(j) only has indeterminate variables that p(i) has, too
Lz = (Ldeg_e * Lz_deg' & ~(~Lvar * Lz_var'));

% discard monomials based on simple checks
[~,Ldegmat] = get_degmat(S);
lz = numel(z);
% perform checks vector-wise
MX = arrayfun(@(i) repmat(max(S.degmat(Ldegmat(i,:),Iv),[],1),lz,1), idx, 'UniformOutput', false);
MN = arrayfun(@(i) repmat(min(S.degmat(Ldegmat(i,:),Iv),[],1),lz,1), idx, 'UniformOutput', false);
% Imx = any(  ceil(max(p.degmat,[],1)/2) < z.degmat , 2 );
% Imn = any( floor(min(p.degmat,[],1)/2) > z.degmat , 2 );
zdm = repmat(z.degmat,lp,1);
Irem = [ceil(vertcat(MX{:})/2) < zdm, floor(vertcat(MN{:})/2) > zdm];
% Lz(:,Imx | Imn) = false;
Lz(reshape(any(Irem,2),lz,lp)') = false;

% remove unused monomials from base vector
I = any(Lz,1);
degmat = z.degmat(I,:);
Lz(:,~I) = [];

% removes monomials outside half Newton polytope
if ~isempty(newton_solver)
    Lz_red = arrayfun(@(i) newton_reduce(S.degmat(Ldegmat(i,:),Iv),degmat,newton_solver), idx, 'UniformOutput', false);
    Lz_red = horzcat(Lz_red{:})';
    
    [Z,K,Mp,Md] = gram_internal(Lz,degmat,z.indets, Lz_red);
else
    [Z,K,Mp,Md] = gram_internal(Lz,degmat,z.indets);	
end

z = (build_monomials(degmat,z.indets)); 

end
