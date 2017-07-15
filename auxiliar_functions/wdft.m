function [Y f] = wdft(X,a,K,fs,Q,D)
% Y = wdft(X,a) encontra a Warped DFT do sinal 'X',onde 'a' é um parâmetro
% (|a| < 1), que indica se a resolução das frequências baixas será maior
% (se a < 0) ou se a resolução das frequências altas será maior (se a > 0).
%
% Y = wdft(X,a,N), se N >length(X), encontra a Overcomplete WDFT (OWDFT) do
% sinal X, ou seja, Y terá tamanho N, maior do que X. Isto é similar à DFT
% com mais bins do que o tamanho original do sinal.
%
% [Y f] = wdft(X,a,[],fs) também retorna o vetor de frequências (warped),
% onde 'fs' indica a frequência de amostragem. Se 'fs' não for
% especificado, é considerado como 1, ou seja, o vetor de frequências 'f' é
% retornado em função de 'fs'.
%
% SÓ FUNCIONA COM LENGTH(X) PAR !

% Última atualização: 23/07/2010

%% Processando as entradas e saídas
is_column = size(X, 1) - 1;
X = X(:);

N = length(X); % Normalmente K é igual a N. Se for menor, é OWDFT

if nargin < 6
    calc_D = true;
    if nargin < 5
        calc_Q = true;
        if nargin < 4
            fs = 1;
            if nargin < 3
                K = N;
                if nargin < 2
                    a = 0;
                end
            end
        end
    else
        calc_Q = false;
    end
else
    calc_Q = false;
    calc_D = false;
end

if isempty(K)
    K = N;
end

if K < N
    error('WDFT - O número de bins deve ser igual ou maior do que o tamanho do vetor de entrada.')
end

if nargout > 1
    f_ = (0:K-1)/K * 2*pi;
    f = 2*atan(((1+a)/(1-a)*tan(f_/2)))*fs/(2*pi);
end

%% Encontrando a matriz Q
% Esta parte pode ser excluída numa eventual implementação em C++, pois
% todas as possíveis matrizes podem estar embutidas no código

if calc_Q
    % Os polinômios são A(z) = -a + z^-1 e Ã(z) = 1 - a*z^-1, ou seja,
    % os coeficientes de A(z) são -a(*z^0) e 1(*z^-1), e   
    % os coeficientes de Ã(z) são 1(*z^0) e -a(*z^-1), e   
    A = cell(1, N);
    A_tilde = cell(1, N);

    % Polinômios elevados a 0 e à 1ª potência
    A{1} = 1;           A{2} = [-a; 1];
    A_tilde{1} = 1;     A_tilde{2} = [1; -a];

    A_1 = A{2};     A_tilde_1 = A_tilde{2}; % Pré-alocação

    % Encontrando os polinômios elevados ao quadrado e acima
    for i = 3:N
        A{i} = conv(A{i-1}, A_1);
        A_tilde{i} = conv(A_tilde{i-1}, A_tilde_1);
    end

    % Encontrando a matriz Qe, coluna a coluna. Repare que Qe depende do
    % tamanho N do vetor de entrada X, e não do número K de pontos da WDFT
    Q = zeros(N); % Q = Qe no caso do mapeamento ser de ordem 1. Se for OWDFT, Q = V*Qe, onde V é uma matriz identidade em cima e com zeros embaixo
    for i = 1:N
        Q(:, i) = conv(A_tilde{N-i+1}, A{i}); % Ã^(N-1-i)*A^i , modificado para que i tenha índice 1 em vez de índice 0
    end
end

%% Computando a WDFT
% Encontrando o vetor P. Se for OWDFT, basta fazer 'zeropad' em vez de
% multiplicar por uma matriz com zeros embaixo
P = fft(Q*X, K);

% Encontrando a matriz diagonal D = diag( 1 ./ fft(Ã^(N-1)) )
if calc_D
    D = fft(A_tilde{N}, K); % Em vez de usar uma matriz diagonal, podemos utilizar a função ./ do Matlab para aplicar elemento a elemento
end

Y = P./D;

if ~is_column
    Y = Y.';
end
