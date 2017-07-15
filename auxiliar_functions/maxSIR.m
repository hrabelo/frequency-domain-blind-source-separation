function [tr perm] = maxSIR(W, Q)
%
% [tr perm] = maxSIR(W, Q) utiliza o método de
% Makino para resolver o problema da permutação em BSS.
% 
% A função encontra o traço de <|W*Q(1:num_senosr)|.^2>, sendo que a ordem das linhas é
% variada, e é retornado o maior valor, assim como a permutação utilizada
% para encontrá-lo.
%
% tr - máximo traço
% perm - vetor coluna com a permutação que obteve o máximo traço
% W - matriz de tamanho N x num_sensors, de mistura
% Q - matriz de tamanho N*num_sensors x num_amostras, com um sinal em cada
%     linha
%
n_src = size(W, 1);
num_sensors = size(W, 2);
N = size(Q,2);

poss_perms = perms(1:n_src);
n_perms = size(poss_perms, 1);
tr = zeros(n_perms, 1);

for prm = 1:n_perms
    fun = zeros(n_src);
    for src = 1:n_src
        b_ind = (src-1)*num_sensors + 1;
        tmp = W(poss_perms(prm, :), :) * Q(b_ind:b_ind + num_sensors - 1, :);
        fun(:, src) = sum(tmp.*conj(tmp), 2) / N;
    end
    tr(prm) = trace(fun);
end
    
[tr ind]= max(tr);
perm = poss_perms(ind, :).';