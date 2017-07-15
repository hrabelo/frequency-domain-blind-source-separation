function indArray = getInvCustomInd(f, dist, f0, fmax)
% indArray = getInvCustomInd(f, dist, f0, fmax) retorna os índices inversos
% da função getCustomInd.
%
% Em outras palavras, para uma frequência f, ela retorna de quais outras
% frequências ela é adjacente ou harmônica.
%
if fmax > f0
    indArray = [getInvCustomHarmInd(f, f0, fmax) getAdjInd(f, dist, f0, fmax)];    % Frequências crescentes
else
    indArray = [getInvCustomHarmInd(f, f0, fmax) getAdjInd(f, dist, fmax, f0)];    % Frequências decrescentes
end     

indArray = unique(indArray);