function indArray = getInvCustomHarmInd(f, f0, fmax)
% indArray = getInvCustomHarmInd(f, f0, fmax) retorna os índices dos harmônicos
% inversos de getCustomHarmInd
%
% Em outras palavras, para uma frequência f, ela retorna de quais outras
% frequências ela é harmônica.

% NÃO ESTÁ FUNCIONANDO MUITO BEM PARA FREQUÊNCIAS DECRESCENTES (pega mais 
% do que o necessário nos ceil, normalmente isso não é problema). Talvez
% tenha que usar floor(.)

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

if f == 1
    indArray = [];
else

    indArray = [2*(f - f0 + 1) + f0,  2*(f - f0 + 1) + f0 - 1, ...
                2*(f - f0) + f0,  2*(f - f0) + f0 - 1, ...
                2*(f - f0 - 1) + f0,  2*(f - f0 - 1) + f0 - 1, ...
                ceil((f - f0 + 1)/2) + f0, ceil((f - f0)/2) + f0, ...
                ceil((f - f0 - 1)/2) + f0];

            if increasing_freq
                indArray = indArray(indArray > f0 & indArray <= fmax);
            else 
                indArray = indArray(indArray >= fmax & indArray < f0);
            end
            
    indArray = unique(indArray);
    indArray(indArray == f) = [];
end

