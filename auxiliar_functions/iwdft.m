function Y = iwdft(X,a,N,invQ,invD)
% Y = iwdft(X,a) encontra a Inverse WDFT do sinal 'X',onde 'a' é o parâmetro
% utilizado na wdft.
%
% Y = iwdft(X,a,N), se N <length(X), encontra a Overcomplete IWDFT (OIWDFT)
% do sinal X, ou seja, Y terá tamanho N, menor do que X.Isto é similar à IDFT
% com mais bins do que o tamanho original do sinal.
%
% SÓ FUNCIONA COM LENGTH(X) PAR !

% Última atualização: 23/07/2010

%% Processando as entradas e saídas
is_column = size(X, 1) - 1;
X = X(:);

K = length(X); % Normalmente N é igual a K. Se for menor, é OIWDFT

if nargin < 5
    calc_invd = true;
    if nargin < 4
        calc_invq = true;
        if nargin < 3
            N = K;
            if nargin < 2
                a = 0;
            end
        end
    else
        calc_invq = false;
    end
else
    calc_invq = false;
    calc_invd = false;
end

if N > K
    error('IWDFT - O número de bins deve ser igual ou maior do que o tamanho do vetor de saída [length(X) >= N].')
end

%% Encontrando a matriz Q^-1
% Esta parte pode ser excluída numa eventual implementação em C++, pois
% todas as possíveis matrizes podem estar embutidas no código

if calc_invq
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

    % Encontrando a matriz Qe, coluna a coluna
    Q = zeros(N); % Q = [Qe ; zeros()] no caso do mapeamento ser de ordem 1, e a IWDFT ser Overcomplete
    for i = 1:N
        Q(:, i) = conv(A_tilde{N-i+1}, A{i}); % Ã^(N-1-i)*A^i , modificado para que i tenha índice 1 em vez de índice 0
    end

    invQ = pinv(Q); % A pseudoinverse de [Qe ; zeros()] é [pinv(Qe) zeros()]. Isto significa ignorar os últimos samples de ifft(invD*X) antes de multiplicar por invQ 
end
    
%% Computando a IWDFT
% Encontrando a matriz diagonal D^-1 = inv(diag( 1 ./ fft(Ã^(N-1)) )) =
% diag( fft(Ã^(N-1) ), sabendo que a inversa de uma matriz diagonal D é 
% igual a 1./D.
if calc_invd
    invD = fft(A_tilde{N}, K); % Podemos utilizar a função .* do Matlab em vez de uma matriz diagonal
end

% x = Q^-1 * transfX, onde transfX = W^-1*D^-1*X
transfX = ifft(invD.*X);
transfX = transfX(1:N, :); % Ignora os últimos samples, para que transfX tenha tamanho N (assim como o da saída)

Y = invQ*transfX;

if ~is_column
    Y = Y.';
end
