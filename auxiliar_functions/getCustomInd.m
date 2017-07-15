function indArray = getCustomInd(f, dist, f0, fmax)
% indArray = getCustomInd(f, dist, f0, fmax) retorna os índices dos harmônicos
% da frequência cujo índice é 'f', e os índices adjacentes a ela, com uma
% distância de 'dist'. 'f0' é o índice da frequência 0, e fmax
% é o índice da frequência máxima.
%
% As frequências podem ser decrescentes. Neste caso, f0 > fmax, ou seja, o
% índice da frequência 0 pode ser 256, e o da frequência máxima 1. Isto é
% útil para encontrar harmônicos de frequências negativas.
%
% A saída é um vetor LINHA
%
% Atualmente ele retorna os seguintes harmônicos:
% f/2 - 1, f/2, f/2 + 1
% 2*f - 1, 2*f, 2*f + 1
%
if fmax > f0
    indArray = [getCustomHarmInd(f, f0, fmax) getAdjInd(f, dist, f0, fmax)];    % Frequências crescentes
else
    indArray = [getCustomHarmInd(f, f0, fmax) getAdjInd(f, dist, fmax, f0)];    % Frequências decrescentes
end     

indArray = unique(indArray);
