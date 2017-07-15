function indArray = getCustomHarmInd(f, f0, fmax)
% indArray = getCustomHarmInd(f, f0, fmax) retorna os índices dos harmônicos
% da frequência cujo índice é 'f', e 'f0' é o índice da frequência 0. fmax
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
if fmax > f0,       increasing_freq = true;             % Frequências crescentes
else                increasing_freq = false;    end     % Frequências decrescentes

if increasing_freq
    if (f < f0) || (f > fmax)
        error('GETCUSTOMHARMIND - Frequência f fora dos limites f0 e fmax')
    end
else
    if (f > f0) || (f < fmax)
        error('GETCUSTOMHARMIND - Frequência f fora dos limites f0 e fmax')
    end
end

if f == f0
    indArray = [];
else

    indArray = [ceil((f-f0)/2) + f0 - 1, ceil((f-f0)/2) + f0, ...
                ceil((f-f0)/2) + f0 + 1, ...                               % f/2 - 1, f/2, f/2 + 1
                (f-f0)*2 + f0 - 1, (f-f0)*2 + f0, (f-f0)*2 + f0 + 1];       % 2*f - 1, 2*f, 2*f + 1

            if increasing_freq
                indArray = indArray(indArray > f0 & indArray <= fmax);
            else 
                indArray = indArray(indArray >= fmax & indArray < f0);
            end
            
    indArray = unique(indArray);
    indArray(indArray == f) = [];
end

