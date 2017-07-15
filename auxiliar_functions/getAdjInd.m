function indArray = getAdjInd(baseind, distance, minind, maxind)
% indArray = getAdjInd(baseind, distance, min, max) retorna os índices adjacentes a
% 'baseind', com uma distância de 'distance', sabendo que o índice mínimo não
% pode ser menor que 'minind' e o índice máximo, maior que 'maxind'
%
% A saída é um vetor LINHA
%
% Exemplos:
% getAdjInd(5, 2, -5, 15) retorna o vetor     [3 4 6 7]
% getAdjInd(5, 4, -5, 15) retorna o vetor [1 2 3 4 6 7 8 9]
% getAdjInd(5, 4, 3, 8) retorna o vetor       [3 4 6 7 8]

indArray = [ (baseind - distance):(baseind - 1) (baseind + 1):(baseind + distance) ];
indArray = indArray(indArray >= minind & indArray <= maxind);