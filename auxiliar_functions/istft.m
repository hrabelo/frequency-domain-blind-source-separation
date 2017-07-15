function x = istft(X, J, varargin)
% x = istft(X, J) 
% x = istft(X, J, maxlen) 
% x = istft(X, J, wind) 
% x = istft(X, J, wind, maxlen) 
% x = istft(X, J, wind, wdft_par) 
% x = istft(X, J, wind, wdft_par, maxlen) 
%         gera a STFT inversa de um ou mais sinais no domínio da frequência,
%         onde
%
%         X - matriz com 3 dimensões, onde os frames são as colunas, cada linha
%             representa um bin de frequência, e na outra dimensão estão 
%             representados os sinais. Estes sinais devem ser SIMÉTRICOS, e
%             se não forem, a função trata como se fossem. 
%             Ver help da stft para mais detalhes
%
%         J - salto entre cada frame da STFT
%
%         maxlen - comprimento máximo da saída. Isto é útil quando se
%                  utilizou "zeropad" na stft, para que a stft inversa seja
%                  igual à entrada da stft. As amostras além do comprimento
%                  máximo são descartadas. 
%                  O default é Inf.
%
%         wind - janela que será utilizada, se aplicável*. Se o tamanho da
%                janela for menor do que o número de bins, isto significa
%                que houve "oversampling"
%
%                * Se foi utilizada uma janela com a propriedade 
%                CONSTANT-OVERLAP-ADD para gerar a STFT (uma Hanning com 50%
%                de overlap, por exemplo), uma janela retangular é suficiente,
%                ou seja, este parâmetro é desnecessário. Por outro lado,
%                alguns codificadores de áudio, por exemplo, utilizam uma
%                janela na STFT e outra na ISTFT, para conseguir a propriedade
%                de CONSTANT-OVERLAP-ADD, que são chamadas de janela de
%                análise e janela de síntese. Neste caso, este parâmetro é
%                essencial.
%
%         x - saída, onde cada linha representa um sinal, e cada coluna um
%             sample

%% Processando as entradas e saídas

% Argumentos padrão
maxlen = Inf;
wdft_par = 0;
wind = ones(1, size(X, 2));
    
if numel(varargin)
    optargin = size(varargin, 2);
    
    switch optargin
        case 0
        case 1
            if length( varargin{1} ) > 1,   wind = varargin{1};  % Se for um vetor, é a janela
            else                            maxlen = varargin{1};   end
        case 2
            wind = varargin{1};
            if varargin{2} >= 1,            maxlen = varargin{2}; % Se não for decimal, é o tamanho máximo
            else                            wdft_par = varargin{2};   end
        case 3    
            wind = varargin{1};
            wdft_par = varargin{2};
            maxlen = varargin{3};
        otherwise
            error('Excesso de argumentos de entrada')
    end
end

N = size(wind, 2); % Tamanho do frame
[M, K, n_frames] = size(X);

if J > N
    error('ISTFT - O salto J não pode ser maior que o tamanho do frame.')
end

if N > K
    error('ISTFT - O tamanho do frame(dado pela janela) não pode ser maior do que o número de bins da DFT (dado pela matriz de entrada X). Diminua o tamanho da janela utilizada.')
end

if N < K
    warning('ISTFT - Oversampling: o tamanho da saída da IDFT é menor que o número de bins de frequência.')
end

%% Reconstruindo o sinal frame a frame
x = zeros(M, (n_frames-1)*J + N);

if wdft_par
    
    A = cell(1, N);
    A_tilde = cell(1, N);
    A{1} = 1;           A{2} = [-wdft_par; 1];
    A_tilde{1} = 1;     A_tilde{2} = [1; -wdft_par];
    A_1 = A{2};     A_tilde_1 = A_tilde{2}; % Pré-alocação
    for i = 3:N
        A{i} = conv(A{i-1}, A_1);
        A_tilde{i} = conv(A_tilde{i-1}, A_tilde_1);
    end
    Q = zeros(N); % Q = [Qe ; zeros()] no caso do mapeamento ser de ordem 1, e a WDFT ser Overcomplete
    for i = 1:N
        Q(:, i) = conv(A_tilde{N-i+1}, A{i}); % Ã^(N-1-i)*A^i , modificado para que i tenha índice 1 em vez de índice 0
    end
    invQ = pinv(Q); % A pseudoinverse de [Q ; zeros()] é [pinv(Q) zeros()]. Isto significa ignorar os últimos samples de ifft(invD*X) antes de multiplicar por invQ 

    invD = fft(A_tilde{N}, K);
    
    for ci = 1:M   
        for frame = 1:n_frames
            x(ci,1 + (frame-1)*J : (frame-1)*J + N) = x(ci,1 + (frame-1)*J : (frame-1)*J + N) + iwdft(X(ci, :, frame), wdft_par, N, invQ, invD) .* wind;   
        end
    end
else
    for ci = 1:M   
        for frame = 1:n_frames
            tmp = ifft(X(ci, :, frame), 'symmetric');
            x(ci,1 + (frame-1)*J : (frame-1)*J + N) = x(ci,1 + (frame-1)*J : (frame-1)*J + N) + tmp(1:N) .* wind;
        end
    end
end

%% Truncando o sinal (se foi utilizado zeropadding na stft)
if ~isinf(maxlen)
    x = x(:, 1:maxlen);
end
