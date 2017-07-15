function [Rf perm] = maxcorr(base_sig, adj_sig)
%
% [Rf, perm] = maxcorr(base_sig, adj_sig) analisa as correlações
% entre um vetor com fontes de uma frequência base e outros vetores com fontes
% de frequências adjacentes, utilizando todas as permutações possíveis das
% fontes da frequência base. É feita a soma dos valores de correlação para cada
% permutação e retornada a maior soma e a permutação que obteve a maior soma.
%
% Rf - máxima soma de correlações
% perm - vetor coluna com a permutação que obteve a máxima soma
% base_sig - matriz com as fontes da frequência base, onde cada fonte é uma linha
% adj_sig - matriz de 3 DIMENSÕES, similar à base_sig, e a 3ª dimensão é a das
%           frequências. É OBRIGATÓRIO que a matriz tenha 3 dimensões
%
N = size(base_sig, 2);
poss_perms = perms(1:size(base_sig, 1)); % Cada linha é uma permutação
n_perms = size(poss_perms, 1);
Rf = zeros(n_perms, 1);

for prm = 1:n_perms
    permbase = base_sig(poss_perms(prm, :), :);
    
    % Pré-calculando para agilizar a função
    med_base = mean(permbase, 2);
    var_base = var(permbase, 0, 2);

    % Pode-se tentar usar a função corrcoef do MatLab, se ela for mais
    % rápida
    for f = 1:size(adj_sig, 1)
        tmpadj = squeeze( adj_sig(f, :, :) );
        Rf(prm) = Rf(prm) + sum( (sum(permbase .* tmpadj, 2)/(N-1) - (N/(N-1)) * med_base .* mean(tmpadj, 2)) ./ (sqrt(var_base .* var(tmpadj, 0, 2))) );
    end
end

[Rf ind]= max(Rf);
perm = poss_perms(ind, :).';
