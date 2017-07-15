function U = dirpattern(W, d, f, teta)
%
% U = dirpattern(W, d, f, teta) retorna o "directivity pattern" de várias
% fontes, onde W é a matriz de separação da fonte. O ângulo perpendicular à
% posição dos sensores é pi/2 (90º)
%
% U - matriz onde cada linha é o "directivity pattern" de uma fonte
%
% teta - vetor linha onde cada elemento é um ângulo onde f deve ser
%        encontrado. Deve estar entre 0 e pi
% W - matriz N (fontes) x M (misturas) separadora
% d - vetor com a posição dos sensores em metros(exemplo [-0.02 0.02], se a
%     distância entre os sensores for 0.04
% f - frequência do sinal
%
global SPEED_OF_SOUND
if isempty(SPEED_OF_SOUND)
    SPEED_OF_SOUND = 343;
end

d = d(:);

c = 1/SPEED_OF_SOUND; % inverso da velocidade de propagação no meio, em s/m
N = size(W, 1); % Número de fontes

U = zeros(N, length(teta));
for cn = 1:N
    U(cn, :) = abs(sum(diag(W(cn, :)) * exp(i*2*pi*f*c*d*cos(teta)), 1));
end