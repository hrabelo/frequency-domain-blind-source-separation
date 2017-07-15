function [X, n_frames] = stft(x, K, J, varargin)
% [X, n_frames] = stft(x, K, J) 
% [X, n_frames] = stft(x, K, J, 'zeropad')
% [X, n_frames] = stft(x, K, J, wind)
% [X, n_frames] = stft(x, K, J, wind, 'zeropad')
% [X, n_frames] = stft(x, K, J, wind, wdft_par) 
% [X, n_frames] = stft(x, K, J, wind, wdft_par, 'zeropad') 
%         gera os frames da STFT de um ou mais sinais no domínio do tempo, 
%         onde
%
%         x - sinal(is) de entrada(s), no domínio do tempo. O sinal deve ser
%             vetores linha, onde as colunas são as amostras, e cada linha é
%             um sinal
%         K - número de 'bins' na frequência. Se não for especificada nenhuma
%             janela, corresponde ao tamanho do frame da STFT 
%         J - salto entre cada frame da STFT
%         wind - janela que será utilizada. É um vetor linha, que será 
%                multiplicado por cada frame antes de aplicar a DFT. Em
%                condições normais, o tamanho de wind é igual a K, a não
%                ser que se deseje "oversampling"
%         'zeropad' - indica que será feito zero-padding. Por default as
%                   últimas amostras serão descartadas (se o número de
%                   amostras do final não for suficiente para gerar um frame
%                   completo). Com este modificador, a função acrescenta
%                   zeros ao final do vetor antes de passar para o domínio
%                   da frequência
%         wdft_par - Se for 0, é aplicada a FFT comum em cada frame, e se
%         for maior que 0, é aplicada a WDFT
%
% Saídas
%     Ao final, é gerada uma matriz X de 3 dimensões, onde os frames são as
%     colunas, cada linha representa um bin de frequência, e na outra dimensão
%     estão representados os sinais
% X =   [ X1(1) X2(1) X3(1) ... Xnum_of_frames(1)
%         X1(2) X2(2) X3(2) ... Xnum_of_frames(2)
%         X1(3) X2(3) X3(3) ... Xnum_of_frames(3)
%                ...                 ...
%         X1(K) X2(K) X3(K) ... Xnum_of_frames(K) ] <- FFT da primeira linha de x
%
%       [ X1(1) X2(1) X3(1) ... Xnum_of_frames(1)
%         X1(2) X2(2) X3(2) ... Xnum_of_frames(2)
%         X1(3) X2(3) X3(3) ... Xnum_of_frames(3)
%                ...                 ...
%         X1(K) X2(K) X3(K) ... Xnum_of_frames(K) ] <- FFT da segunda linha de x
%
% n_frames é o número de frames gerado
%
%
% ATENÇÃO! Algumas amostras do final do vetor x de entrada podem ser
% descartadas no processo, ou seja, a DFT inversa da matriz X NÃO GERARÁ
% exatamente o vetor x. Para evitar o descarte, utilize 'zeropad'.

% Última atualização: 23/07/2010

%% Processando as entradas e saídas

% Argumentos padrão
zeropad = 0;
wdft_par = 0;
wind = ones(1, K);
    
if numel(varargin)
    if find( strcmp('zeropad', varargin) , 1),  zeropad = 1;        end

    optargin = size(varargin, 2) - zeropad; % zeropad é sempre o último argumento, se existir
    
    if optargin > 0
        wind = varargin{1};
        if optargin > 1
            wdft_par = varargin{2};
            if optargin > 2
                error('Excesso de argumentos de entrada')
            end
        end
    end
end

N = size(wind, 2); % Tamanho do frame
M = size(x, 1); % número de sinais
num_samp = size(x, 2); % número de amostras

if J > N
    error('STFT - O salto J não pode ser maior que o tamanho do frame.')
end

if N > K
    error('STFT - O tamanho do frame não pode ser maior do que o número de bins da DFT. Diminua o tamanho da janela utilizada.')
end

if N < K
    warning('STFT - Oversampling: o tamanho da DFT é maior que o número de amostras por frame.')
end

%% Zero-Padding
if zeropad
% Se mod(num_samp-N, J for maior que zero, sobrou um resto, então o número
% de frames deve ser igual a floor((num_samp - N) / J) + 2, e devem ser 
% adicionados zeros ao final da amostra para que a conta feche. O número de
% zeros adicionados é o tamanho do pulo menos o resto, ou seja, 
% J-mod(num_samp - N, J)
    if mod(num_samp - N, J)
        x = [x zeros(M, J - mod(num_samp - N, J))];

% Se o mod(num_samp - N, J) for zero, quer dizer que nenhuma amostra será
% descartada, e o número de frames é igual a floor((size(x,2) - N) / J) + 1,
% ou seja, zeropad deve ser 0
    else
        zeropad = 0;
    end
end

%% Aplicando a FFT ou WDFT frame a frame
n_frames = floor((num_samp - N) / J) + 1 + zeropad;
X = zeros(M, K, n_frames);
if wdft_par
    
    A = cell(1, N);
    A_tilde = cell(1, N);
    A{1} = 1;           A{2} = [-wdft_par; 1];
    A_tilde{1} = 1;     A_tilde{2} = [1; -wdft_par];
    A_1 = A{2};     A_tilde_1 = A_tilde{2};
    for i = 3:N
        A{i} = conv(A{i-1}, A_1);
        A_tilde{i} = conv(A_tilde{i-1}, A_tilde_1);
    end
    Q = zeros(N);
    for i = 1:N
        Q(:, i) = conv(A_tilde{N-i+1}, A{i});
    end
    
    D = fft(A_tilde{N}, K);
    
    for ci = 1:M
        for frame = 1:n_frames
            X(ci, :, frame) = wdft(x( ci, 1 + (frame-1)*J : (frame-1)*J + N ) .* wind, wdft_par, K, 1, Q, D).';
        end 
    end
else
    for ci = 1:M
        for frame = 1:n_frames
            X(ci, :, frame) = fft(x( ci, 1 + (frame-1)*J : (frame-1)*J + N ) .* wind, K).';
        end 
    end
end
