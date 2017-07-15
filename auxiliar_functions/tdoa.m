function r = tdoa(W, ref, f)
%
% r = tdoa(W, ref, f) retorna o Time Difference of Arrival
%
% r - matriz N x (M-1), onde cada linha é o TDOA de uma fonte em relação
%     aos microfones
% W - matriz N (fontes) x M (misturas) separadora
% ref - microfone a ser utilizado como referência. Ref é o número da coluna
%       de W que corresponde a este sensor
% f - frequência do sinal
%
N = size(W,1);
M = size(W,2);

% Exceção
if f == 0
    r = nan(N,M-1);
    return
end

if N == M,      A = inv(W);
else            A = pinv(W);                                        end

r = zeros(N,M-1);
for cm = 1:ref-1
    r(:,cm) = (angle(A(cm,:) ./ A(ref,:)) / (2*pi*f)).'; % Pode-se colocar um sinal de '-' aqui, o que altera apenas a referência
end
for cm = ref+1:M
    r(:,cm-1) = (angle(A(cm,:) ./ A(ref,:)) / (2*pi*f)).'; % Pode-se colocar um sinal de '-' aqui, o que altera apenas a referência
end