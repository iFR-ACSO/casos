function [Ksdp_g_s, Qgram_g, Zgram_g] = zero_diagonal(g, Qgram_g, Zgram_g, Ksdp_g_s, Ml)
% zero diagonal algorithm (similar in practive to simplify from sosopt)

while true 

% replace sum-of-squares constraints by equality
gdiff = (g - [zeros(Ml,1); casos.PS(Zgram_g,Qgram_g)]);

% handle (new) equality constraints
[Qdiff_g,~] = poly2basis(gdiff);

% index of diagonal terms in Q
idx = linspace(1,Ksdp_g_s.^2, Ksdp_g_s);
Qdiag = Qgram_g(idx);

% SCATCHY: use jacobian to verify where does each diagonal term appear 
% in the constraints 
MJ = casadi.DM(jacobian(Qdiff_g,Qdiag));
MJ = MJ.full();
idx = arrayfun(@(i) find(MJ(:,i)), 1:Ksdp_g_s, 'UniformOutput', false);

flag = 0;
for i=1:length(idx)
    if length(symvar(Qdiff_g(idx{i})))==1
        fun = casadi.Function('fun', {Qdiag(i)}, {Qdiff_g(idx{i})});
        if (full(fun(0)) == 0)
            % flag to identify that a change has been made
            flag = 1;

            % I need to remove column and row
            Qgram_g = Qgram_g.reshape(Ksdp_g_s,Ksdp_g_s);
            Zgram_g = Zgram_g.reshape(Ksdp_g_s,Ksdp_g_s);

            % i is the column and row
            columnsToKeep = [1:(i-1), (i+1):length(idx)];
            
            % update Qgram_g
            Qgram_g = Qgram_g(:,columnsToKeep);
            Qgram_g = Qgram_g(columnsToKeep, :);
            Qgram_g = Qgram_g(:);

            % update Zgram_g
            Zgram_g = Zgram_g(:, ~ismember(1:size(Zgram_g, 2), i));
            Zgram_g = Zgram_g(~ismember(1:size(Zgram_g, 2), i), :);
            Zgram_g = Zgram_g(:);
            Ksdp_g_s = Ksdp_g_s-1;
            break
        end
    end
end

% exit if not changes took place in the for
if flag==0, break; end

end
end