function teta = doa(W, d, f)
%
% teta = doa(W, d, f) retorna o ângulo das fontes ou NaN se o ângulo não 
% puder ser encontrado (o resultado do arccos é um número complexo)
%
% teta - vetor coluna onde cada elemento é o ângulo de uma fonte
% W - matriz N (fontes) x M (misturas) separadora
% d - espaçamento entre os sensores, em metros
% f - frequência do sinal
%
% OS SENSORES DEVEM ESTAR DISPOSTOS EM UMA LINHA E IGUALMENTE ESPAÇADOS
%
global SPEED_OF_SOUND
if isempty(SPEED_OF_SOUND)
    SPEED_OF_SOUND = 343;
end

c = 1/SPEED_OF_SOUND; % inverso da velocidade de propagação no meio, em s/m
N = size(W,1);

% Exceção
if f == 0
    teta = nan(N,1);
    return
end

combs_2by2 = nchoosek(1:size(W,2), 2);
n_cmb = size(combs_2by2, 1);

% Nota: a expressão (W^-1)*(P^-1) = (W^-1)*P', como P é uma matriz de
% permutação, apenas altera as colunas de W^-1. Como existe um ângulo para
% cada coluna(cada coluna de W^-1 representa uma fonte, diferente de W,
% onde cada coluna representa uma mistura), a permutação não altera o
% resultado final
A = inv(W); % A = W^-1

% Testa para todas as combinações, para o caso de alguma divisão der maior
% do que 1 em módulo
for cmb = 1:n_cmb
    mic1 = combs_2by2(cmb, 1);
    mic2 = combs_2by2(cmb, 2);
    tmp = angle(A(mic1,:) ./ A(mic2,:)) / (2*pi*f*c*(mic2 - mic1)*d); % O (mic2 - mic1) gera a distância d correta
    
    % Como acos de um número maior que 1 dá um número complexo, se nenhum
    % dos valores for maior do que 1, pode sair do loop
    if ~sum(abs(tmp) > 1)
        break;
    end
end
    
if ~sum(abs(tmp) > 1)
    teta = acos(tmp).';
else
    teta = nan(N, 1);
end