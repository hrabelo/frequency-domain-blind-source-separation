function [P ind_solve] = bss_fdpermsolve(W, Y, varargin)
%
%   P = bss_fdpermsolve(W, Y)
%
%   W - Matriz com 3 dimensões, onde a 3ª dimensão é a das frequências, e
%   as linhas e colunas são os coeficientes da matriz de desmistura. Cada
%   linha corresponde aos coeficientes das misturas que gerarão uma fonte
%
%   Y - Matriz com as saídas, ou seja, W*X, onde X é a matriz das misturas,
%   em cada frequência. Cada linha é uma fonte, cada coluna é um frame,
%   e a 3ª dimensão é a das frequências
%
%   O tamanho da FFT deve ser par!!!
%
%       Properties:
%               'Method' - modifica o método de resolver o problema da
%               permutação. Pode ser:
%                   'tdoa'
%                   'conjcorr' (default)
%                   'harmcorr'
%                   'globalcorr'
%                   'doa'
%                   'doa_adjcorr'
%                   'doa_harmcorr'
%                   'doa_conjcorr'
%                   'supervised' - necessita do sinal da fonte
%                   'maxsir' - necessita do sinal da fonte em cada sensor
%
%               'useSymmetry' - explora a simetria da FFT, se existir. Por
%               default, não explora a simetria.
%
%           USANDO DOA
%               'DoaThreshold' - diferença máxima entre o ângulo 
%               encontrado no DOA para uma frequência e a média dos ângulos
%               para determinar se a medida é confiável. O padrão é
%               1.5*std(teta).
%
%               'DirPatternThreshold' - diferença mínima entre os
%               'directivity patterns' dos ângulos das fontes, em dB.
%
%               'DistanceBetweenSensors' - distância entre os microfones
%               para utilizar no método DOA, em metros. OBRIGATÓRIO
%
%               'SampFrequency' - frequência de amostragem.
%
%           USANDO TDOA
%
%               'useSymmetry' é OBRIGATÓRIO
%
%               'MaxDistBetweenSensors' - distância máxima entre os microfones
%               para utilizar no método TDOA, em metros. Utilizada apenas
%               para limitar as frequências onde a técnica é aplicada. Se
%               for omitida, não há limite.
%
%           USANDO CORRELAÇÃO
%               'NumAdjFrequencies' - número de frequências adjacentes
%               utilizadas no método da correlação, para um lado e para o
%               outro. A distância (em Hz) é igual a
%                   NumAdjFrequencies*fs/Nfreqs, onde Nfreqs=size(W, 1)
%               e o número total de frequências é igual a 
%                   2*NumAdjFrequencies
%               O padrão é 3
%
%               'CorrThreshold' - valor mínimo de correlação para que o
%               método de correlação seja confiável. Só é aplicável se
%               também estiver sendo utilizado o método de correlação
%               harmônica. O padrão é 
%                   Nfontes*NumAdjFrequencies, onde Nfontes=size(W, 2)
%
%               'HarmonicThreshold' - valor mínimo de correlação para que o
%               método de correlação harmônica seja confiável. O padrão é
%                   Nfontes*NumHarmFrequencies*0,5
%
%               'Envelope' - tipo de envelope, 'AbsValue', 'PowValue',
%                            'PowValue2'
%
%           USANDO MODO SUPERVISIONADO
%               'SourceSignal' - sinal original da fonte, no mesmo formato
%               de Y, ou seja, a dimensão 3 é a das frequências, e cada
%               linha é o sinal de uma fonte
%
%           USANDO MODO MAXSIR
%               'MicSourceComponents' - sinal com os sinais originais das
%               fontes, assim como vistos em cada microfone, ou seja, após
%               passar pela função de transferência da sala. O formato é o
%               mesmo de Y, mas a dimensão 2 (as linhas) é distribuída da 
%               seguinte forma: 
%               linha(1)            Source1_in_Mic1 
%               linha(2)            Source1_in_Mic2
%                                       ...
%               linha(M)            Source1_in_MicM
%               linha(M+1)          Source2_in_Mic1 
%               linha(M+2)          Source2_in_Mic2
%                                       ...
%               linha(2*M)          Source2_in_MicM
%                                       ...
%               linha[(N-1)*M+1]    SourceN_in_Mic1 
%               linha[(N-1)*M+2]    SourceN_in_Mic2
%                                       ...
%               linha(N*M)          SourceN_in_MicM

%% Inicialização
debug = 1;

global SPEED_OF_SOUND
SPEED_OF_SOUND = 343;
KMEANS_FACTOR = 0.5;
KMEANS_MAXITER = 5;
tdoa_cluster = 'custom';

K = size(W, 1); % Número de bins de frequência
N = size(W, 2); % Número de fontes
M = size(W, 3); % Número de misturas
num_of_frames = size(Y, 3);
c_ind_solve = 1; % Índice da saída


%% Processando as entradas e saídas
if numel(varargin)

    ind_arg1 = find( strcmp('Method', varargin) , 1);
    ind_arg2 = find( strcmp('DistanceBetweenSensors', varargin) , 1);
    ind_arg3 = find( strcmp('SourceSignal', varargin) , 1);
    ind_arg4 = find( strcmp('DoaThreshold', varargin) , 1);
    ind_arg5 = find( strcmp('SampFrequency', varargin) , 1);
    ind_arg6 = find( strcmp('NumAdjFrequencies', varargin) , 1);
    ind_arg7 = find( strcmp('CorrThreshold', varargin) , 1);
    ind_arg8 = find( strcmp('HarmonicThreshold', varargin) , 1);
    ind_arg9 = find( strcmp('MicSourceComponents', varargin) , 1);
    ind_arg10 = find( strcmp('useSymmetry', varargin) , 1);
    ind_arg11 = find( strcmp('DirPatternThreshold', varargin) , 1);
    ind_arg12 = find( strcmp('MaxDistBetweenSensors', varargin) , 1);
    ind_arg13 = find( strcmp('Envelope', varargin) , 1);
     
    if ~isempty(ind_arg1),      method = varargin{ind_arg1 + 1};
    else                        method = 'conjcorr';                end
    if ~isempty(ind_arg2),      d = varargin{ind_arg2 + 1};         end
    if ~isempty(ind_arg3),      S = varargin{ind_arg3 + 1};         end
    if ~isempty(ind_arg4),      doa_thre = varargin{ind_arg4 + 1};  
    else                        doa_thre = 0;                       end
    if ~isempty(ind_arg5),      fs = varargin{ind_arg5 + 1};        end
    if ~isempty(ind_arg6),      corr_neighbor_freq = varargin{ind_arg6 + 1};
    else                        corr_neighbor_freq = 3;             end
    if ~isempty(ind_arg7),      corr_thre = varargin{ind_arg7 + 1};
    else                        corr_thre = N * corr_neighbor_freq; end  %%%%% ARRUMAR %%%%%%%%%%%%%%
    if ~isempty(ind_arg8),      harmonic_thre = varargin{ind_arg8 + 1};
    else                        harmonic_thre = 0;                  end  %%%%% ARRUMAR %%%%%%%%%%%%%%
    if ~isempty(ind_arg9),      Q = varargin{ind_arg9 + 1};         end
    if ~isempty(ind_arg10),     unique_K = K/2 + 1; % Número de bins únicos
    else                        unique_K = K;                       end
    if ~isempty(ind_arg11),     u_thre = varargin{ind_arg11 + 1};   
    else                        u_thre = -Inf;                      end
    if ~isempty(ind_arg12),     dmax = varargin{ind_arg12 + 1};   
    else                        dmax = 0;                           end
    if ~isempty(ind_arg13),     env_str = varargin{ind_arg13 + 1};
    else                        env_str = 'AbsValue';               end    

else
    % DEFAULT
    method = 'conjcorr';
    corr_neighbor_freq = 3;
    env_str = 'AbsValue';
end

switch env_str
    case 'AbsValue',        env_type = 1;
    case 'PowValue',        env_type = 2;
    case 'PowValue2',       env_type = 3;
    case 'SPDValue',        env_type = 4;
end


%% Mais inicialização
must_do_doa = false;
doa_no_tests = false;
global_corr_last = false;
must_do_tdoa = false;
must_do_envprecalc = false;
must_do_precalc = false;
must_do_adjprecalc = false;
must_do_corrprecalc = false;
must_do_corr = false;
must_do_localcorr = false;
must_do_harmonic = false;
must_do_precorr = false;
must_do_globalcorr = false;
must_do_supervised = false;
must_do_maxSIR = false;
corr_conf_condition = false;

poss_perms = perms(1:N); % Cada linha é uma permutação possível
n_perms = size(poss_perms, 1); % Número de permutações

P = (1:N).' * ones(1, K); % Matriz de permutação. Cada coluna é uma permutação
Ptmp = (1:N).' * ones(1, unique_K); % Matriz de permutação temporária. Por enquanto, é igual a P nas frequências confiáveis e diferente nas
                                    % não-confiáveis. Futuramente, é melhor que seja válida apenas dentro de cada método.
                                    % MODIFICAR  MODIFICAR  MODIFICAR MODIFICAR  MODIFICAR  MODIFICAR MODIFICAR  MODIFICAR  MODIFICAR MODIFICAR 
ind_conf = logical(false(1, unique_K));
envY = zeros(unique_K, N, num_of_frames);

%% Escolha do método
switch lower(method)
    case 'conjcorr'
        must_do_envprecalc = true;
        must_do_precalc = true;
        must_do_corrprecalc = true;
        must_do_corr = true;
        ind_solve = cell(1,1);

    case 'globalcorr'
        must_do_envprecalc = true;
        must_do_precalc = true;
        must_do_corrprecalc = true;
        must_do_globalcorr = true;
        global_corr_last = true;
        ind_solve = cell(1,1);
        
    case 'localcorr'
        must_do_envprecalc = true;
        must_do_precalc = true;
        must_do_corrprecalc = true;
        must_do_localcorr = true;
        ind_solve = cell(1,1);

    case 'allcorr'
        must_do_envprecalc = true;
        must_do_precalc = true;
        must_do_corrprecalc = true;
        must_do_globalcorr = true;
        must_do_localcorr = true;
        ind_solve = cell(1,1);
        
    case 'doa'
        must_do_doa = true;
        doa_no_tests = true;
        ind_solve = cell(1,1);
        
    case 'doa_adjcorr'
        must_do_doa = true;
        must_do_envprecalc = true;
        must_do_precalc = true;
        must_do_adjprecalc = true;
        must_do_corr = true;
        corr_conf_condition = true;
        ind_solve = cell(2,1);
        
    case 'doa_harmcorr'
        must_do_doa = true;
        must_do_envprecalc = true;
        must_do_precalc = true;
        must_do_adjprecalc = true;
        must_do_precorr = true;
        must_do_harmonic = true;
        must_do_corr = true;
        corr_conf_condition = true;
        ind_solve = cell(4,1);

    case 'harmcorr'
        must_do_envprecalc = true;
        must_do_precalc = true;
        must_do_adjprecalc = true;
        must_do_precorr = true;
        must_do_harmonic = true;
        must_do_corr = true;
        corr_conf_condition = true;
        ind_solve = cell(3,1);
        
    case 'tdoa_conjcorr'
        must_do_tdoa = true;
        must_do_envprecalc = true;
        must_do_precalc = true;
        must_do_corrprecalc = true;
        must_do_corr = true;
        corr_conf_condition = true;
        ind_solve = cell(2,1);
    
    case 'doa_conjcorr'
        must_do_doa = true;
        must_do_envprecalc = true;
        must_do_precalc = true;
        must_do_corrprecalc = true;
        must_do_corr = true;
        corr_conf_condition = true;
        ind_solve = cell(2,1);
 
    case 'doa_globalcorr'
        must_do_doa = true;
        must_do_envprecalc = true;
        must_do_precalc = true;
        must_do_corrprecalc = true;
        must_do_globalcorr = true;
        global_corr_last = true;
        ind_solve = cell(2,1);
        
    case 'doa_allcorr'
        must_do_doa = true;
        must_do_envprecalc = true;
        must_do_precalc = true;
        must_do_corrprecalc = true;
        must_do_globalcorr = true;
        must_do_localcorr = true;
        ind_solve = cell(2,1);
        
    case 'tdoa'
        must_do_tdoa = true;
        ind_solve = cell(1,1);
        
    case 'supervised'
        must_do_supervised = true;
        ind_solve = cell(1,1);
        
    case 'maxsir'
        must_do_maxSIR = true;
        ind_solve = cell(1,1);
end

%% Envelope
if must_do_envprecalc
    if debug
        disp('--- Pré-cálculos do envelope para os métodos de correlação:')
    end
    
    for k = 1:unique_K
        Ytmp(:, :) = Y(k, :, :);
        envY(k, :, :) = envelope(Ytmp, env_type, squeeze(W(k, :, :)));
    end
    
end

%% DOA
if must_do_doa
    if debug
        disp('--- Método DOA (Direction of Arrival):')
    end
    
    freq = [0 : fs/K : (fs/K)*(K/2-1), -fs/2 : fs/K : -fs/K];
    mics = -(M-1)*d/2 : d : (M-1)*d/2;
    teta = zeros(N, unique_K);
    
    for k = find(~ind_conf)
        teta(:, k) = doa( squeeze(W(k, :, :)), d, freq(k) );
        [teta(:, k), Ptmp(:, k)] = sort(teta(:, k), 'ascend');
    end

    % Testa se o ângulo pode ser encontrado
    ind_conf(sum(isnan(teta), 1) == 0) = true;
    mteta = median(teta(:, ind_conf), 2);
    
    if(~doa_thre), doa_thre = std( teta(:, ind_conf), 0, 2) * 1.5;  end
    
    % Testa se o ângulo encontrado é próximo da média dos ângulos
    if ~doa_no_tests
        for k = find(ind_conf)
            if sum( abs(teta(:, k) - mteta) > doa_thre) > 0
                ind_conf(k) = false;
            end
        end
    end
    
    % Testa se o SIR entre os "directivy pattern" dos ângulos das fontes é
    % alto
    if ~doa_no_tests
        for k = find(ind_conf)
            Wtmp = squeeze(W(k, Ptmp(:, k), :)); % Ajusta a permutação de W

            % A função abaixo tem dimensão num_fontes x num_fontes. Cada linha
            % representa o "directivity pattern" em db de uma das fontes para os
            % ângulos estimados para todas as fontes. A diagonal principal
            % representa a amplitude U de uma fonte no ângulo estimado para
            % ela, e os elementos fora da diagonal principal representam a
            % amplitude U desta fonte em outros ângulos, que devem ser
            % maiores, se a estimativa estiver correta (a diagonal
            % principal é o mínimo)
            fun = 10*log10( dirpattern( Wtmp, mics, freq(k), teta(:, k).' ).^2 );

            fun = sum(sum(triu(fun, 1))) + sum(sum(tril(fun, -1))) - trace(fun);
            if(fun < u_thre)
                ind_conf(k) = false;
            end
        end
    end

    P(:, ind_conf) = Ptmp(:, ind_conf);
    
    % Atualizando as matrizes
    for k = find(ind_conf)
        W(k, :, :) = W(k, P(:, k), :);
        envY(k, :, :) = envY(k, P(:, k), :);
    end
    
    ind_solve{c_ind_solve} = ind_conf;
    c_ind_solve = c_ind_solve + 1;
    
    if debug
        disp( sprintf('Número de bins resolvidos: %.0f', sum(ind_conf)) )
        disp('Desvio Padrão dos ângulos estimados (em radianos):')
        disp(doa_thre/1.5)
        disp('Ângulos encontrados (em graus):')
        disp(mteta*180/pi)
    end
end

%% TDOA
if must_do_tdoa
    if debug
        disp('--- Método TDOA (Time Difference of Arrival):')
    end
    
    freq = [0 : fs/K : (fs/K)*(K/2-1), -fs/2 : fs/K : -fs/K];

    % Se a condição abaixo for verdadeira, o TDOA funciona para todas as
    % frequências
    if dmax < SPEED_OF_SOUND/fs,         dmax = 0;                  end

    % Encontra a frequência até a qual o algoritmo funciona
    if dmax,    ind_fmax = find(freq > SPEED_OF_SOUND/(2*dmax), 1) - 1;
    else        ind_fmax = length(ind_conf);                        end
    
    r = zeros(N, ind_fmax, M-1);
    ref_mic = 2;

    for k = 1:ind_fmax
        if ~ind_conf(k)
            r(:, k, :) = tdoa(squeeze(W(k, :, :)), ref_mic, freq(k));
        else
            r(:, k, :) = nan(N, M-1);
        end
    end
    
    ind_valid = false(size(ind_conf));
    ind_valid(~isnan(r(1, :, 1))) = true; % Lógico, 0 se nan, 1 se valor válido
    r = r(:, ind_valid, :); % Limpa a matriz 'r', para não dar Warning no kmeans
    
    map_freq2r = zeros(size(ind_valid));
    map_freq2r(ind_valid) = 1:size(r, 2); % map_freq2r contém zeros onde r era 'nan' Ex.: [0 0 1 2 3 0 4 5 6 7 0 8]
    freq_valid = find(ind_valid); % freq_valid diz os índices das frequências que não contém 'nan'
    switch lower(tdoa_cluster)
        case 'custom'
            perm = zeros(size(Ptmp));
            tmp = zeros(1, n_perms);
            tmp2 = zeros(N, M-1);
            tmpr = r; % Contém r depois de permutado. Serve para encontrar os centróides
            centr = zeros(N, M-1);
            converged = false;

            while(~converged)
                centr(:, :) = mean(tmpr, 2);  % Encontra os N centróides. A dimensão é N x (M-1)
                
                for k = find(~ind_conf & ind_valid)
                    rk =  map_freq2r(k); % A frequência k é correspondente ao índice rk na matriz r 
                    if ~rk,              continue;       end % Se a frequência deu 'nan' no tdoa

                    % Testa as permutações e retorna a soma das distâncias
                    % de cada uma
                    for i = 1:n_perms
                        tmp2(:, :) = r(poss_perms(i, :), rk, :);
                        tmp(i) = sum(  sum((tmp2 - centr).^2, 2)  ); % A soma interna é a norma^2 da diferença entre os vetores
                    end

                    [dummy tmpind] = min(tmp); % A menor distância simboliza a permutação correta
                    perm(:, k) = poss_perms(tmpind, :).';
                end

                % Testando convergência
                remaining_freqs = sum(sum(abs(perm(:, ~ind_conf & ind_valid) - Ptmp(:, ~ind_conf & ind_valid)), 1) ~= 0);
                if ~remaining_freqs
                    converged = true;
                end
                
                for k = find(~ind_conf & ind_valid)
                    Ptmp(:, k) = perm(:, k); % Atualiza Ptmp
                    tmpr(:, rk, :) = r(Ptmp(:, k), rk, :); % Permuta r de acordo com o resultado
                end
            end
            
%            r = permute(r, [1 3 2]);
%            cluster_sum = sum(  sum(bsxfun(@minus, r, centr).^2, 2)  , 3);
%            cluster_dist = squeeze(sum((bsxfun(@minus, r, centr)).^2, 2));
%            tdoa_thre = std( cluster_dist(:, ind_conf), 0, 2) * 1;
%            tdoa_conf = (sum(bsxfun(@gt, cluster_dist, tdoa_thre), 1) == 0);
            ind_conf(~ind_conf & ind_valid) = true;
            
        case 'kmeans'
            r = reshape(r, size(r, 1)*size(r, 2), M-1);
            
            cluster_valid = false;
            cc = 0;
            while ~cluster_valid && cc < KMEANS_MAXITER

                [dummy centr cluster_sum cluster_dist] = kmeans(r, N); % dist é uma matriz com o número de linhas de r e N colunas

                cluster_valid = true;
                for cn = 1:N
                    if sum(dummy == cn) < KMEANS_FACTOR*(length(dummy)/N)
                        cluster_valid = false;
                    end
                end
                cc = cc + 1;
            end

            % Cada centróide simboliza uma fonte real
            tmpperm = zeros(N, 1);
            for cc = 1:length(freq_valid)
                tmpmat = cluster_dist((cc-1)*N + 1 : (cc-1)*N + N, :); % Matriz fontes permutadas x fontes reais

                for cn = 1:N
                    [dummy perm_srcs] = min(tmpmat); % Vetor linha com as fontes mais próximas de cada centróide
                    [dummy real_src]= min(dummy); % Escalar com o índice do centróide com o menor valor, i.e, a fonte real  

                    tmpperm(perm_srcs(real_src)) = real_src; % Atualiza a permutação
                    tmpmat(perm_srcs(real_src), :) = inf(1,N); % Elimina a linha e coluna da matriz
                    tmpmat(:, real_src) = inf(N,1);
                end

                Ptmp(:, freq_valid(cc)) = tmpperm;
                ind_conf(freq_valid(cc)) = true; % MODIFICAR  MODIFICAR  MODIFICAR  MODIFICAR  MODIFICAR  MODIFICAR  MODIFICAR  MODIFICAR  MODIFICAR 
            end
            
    end

    % Atualizando as matrizes
    for k = find(ind_conf)
        P(:, k) = Ptmp(:, k);
        W(k, :, :) = W(k, P(:, k), :);
        envY(k, :, :) = envY(k, P(:, k), :);
    end
    
    ind_solve{c_ind_solve} = ind_conf;
    c_ind_solve = c_ind_solve + 1;
    
    if debug
        disp( sprintf('Número de bins resolvidos: %.0f', sum(ind_conf)) )
%        disp('Soma quadrática de cada cluster:')
%        disp(cluster_sum)
        disp('Centróides encontrados:')
        disp(centr)
    end
end

%% Pré-cálculos para otimização
if must_do_precalc
    if debug
        disp('--- Pré-cálculos para otimização:')
    end

    Ymean = zeros(N, unique_K);
    Yvar = zeros(N, unique_K);
    
    if debug
        disp('------ Médias e Variâncias')
    end

    Ytmp = zeros(N, num_of_frames);
    for k = 1:unique_K
        Ytmp(:, :) = envY(k, :, :);
        Ymean(:, k) = mean(Ytmp, 2);
        Yvar(:, k) = var(Ytmp, 0, 2);
    end
end

%% Pré-cálculos para otimizar o método de correlação adjacente
if must_do_adjprecalc
    if debug
        disp('--- Pré-cálculos para agilizar o método de correlação de frequências adjacentes:')
    end

    Ycorr = repmat(struct( 'indFreqs', [], 'AdjFreqCorr', []), 1, unique_K);

    Ytmp = zeros(N, num_of_frames);
    Ytmp2 = zeros(N, num_of_frames);
    for k = find(~ind_conf)
        Ytmp(:, :) = envY(k, :, :);
        if(k <= K/2 + 1) % Se a frequência for positiva ou -fs/2 ( se 'useSymmetry' estiver habilitado, é necessário fazer isto. Se não estiver, como o envelope não depende da fase, envelope(Y_fs/2) = envelpe(Y_-fs/2) ))
            kk = getAdjInd(k, corr_neighbor_freq, 1, K/2 + 1);
        else% Se a frequência for negativa
            kk = getAdjInd(k, corr_neighbor_freq, K/2 + 1, unique_K);
        end            
        
        Ycorr(k).indFreqs = kk;
        Ycorr(k).AdjFreqCorr = zeros(length(kk), N, N);
        for i = 1:length(kk)
            Ytmp2(:, :) = envY(kk(i), :, :);
            Ycorr(k).AdjFreqCorr(i, :, :) = fast_corr2by2(Ytmp, Ytmp2, Ymean(:, k), Ymean(:, kk(i)), Yvar(:, k), Yvar(:, kk(i)));
        end
    end
    
end

%% Pré-cálculos para otimizar o método de correlação conjugada
if must_do_corrprecalc
    if debug
        disp('--- Pré-cálculos para agilizar o método de correlação de frequências:')
    end
    
    Ycorr = repmat(struct( 'indFreqs', [], 'AdjFreqCorr', []), 1, unique_K);

    Ytmp = zeros(N, num_of_frames);
    Ytmp2 = zeros(N, num_of_frames);
    for k = find(~ind_conf)
        Ytmp(:, :) = envY(k, :, :);
        if(k <= K/2 + 1) % Se a frequência for positiva ou -fs/2 ( se 'useSymmetry' estiver habilitado, é necessário fazer isto. Se não estiver, como o envelope não depende da fase, envelope(Y_fs/2) = envelpe(Y_-fs/2) ))
            kk = getCustomInd(k, corr_neighbor_freq, 1, K/2 + 1);
        else% Se a frequência for negativa
            kk = getCustomInd(k, corr_neighbor_freq, K/2 + 1, unique_K);
        end            
        
        Ycorr(k).indFreqs = kk;
        Ycorr(k).AdjFreqCorr = zeros(length(kk), N, N);
        for i = 1:length(kk)
            Ytmp2(:, :) = envY(kk(i), :, :);
            Ycorr(k).AdjFreqCorr(i, :, :) = fast_corr2by2(Ytmp, Ytmp2, Ymean(:, k), Ymean(:, kk(i)), Yvar(:, k), Yvar(:, kk(i)));
        end
    end
    
end

%% Correlação global
if must_do_globalcorr
    if debug
        disp('--- Correlação global:')
    end

    Ptmp = (1:N).' * ones(1, unique_K); % Inicializa Ptmp. Pesquisar uma forma melhor de fazer isso...
    perm = Ptmp;
    
    centY = zeros(N, num_of_frames);
    tmpenvY = envY; % Serve para calcular a média
    Ytmp = zeros(N, num_of_frames);
    Ytmp2 = zeros(unique_K, num_of_frames);
    tmp = zeros(1, n_perms);
    converged = false;

    while(~converged)
        for cn = 1:N
            Ytmp2(:, :) = tmpenvY(:, cn, :);
            centY(cn, :) = mean(Ytmp2, 1);
        end
        centmean = mean(centY, 2);
        centvar = var(centY, 0, 2);

        for k = find(~ind_conf)
            Ytmp(:, :) = envY(k, :, :);
            tmpRf = fast_corr2by2(Ytmp, centY, Ymean(:, k), centmean, Yvar(:, k), centvar);

            for i = 1:n_perms
                tmp(i) = sum(diag(tmpRf(poss_perms(i, :), :)));
            end

            [dummy tmpind] = max(tmp);
            perm(:, k) = poss_perms(tmpind, :).';
        end

        % Testando convergência
        remaining_freqs = sum(sum(abs(perm(:, ~ind_conf) - Ptmp(:, ~ind_conf)), 1) ~= 0);
        if ~remaining_freqs
            converged = true;
        end

        % Atualizando as matrizes
        for k = find(~ind_conf)
            Ptmp(:, k) = perm(:, k); % Atualiza Ptmp
            tmpenvY(k, :, :) = envY(k, Ptmp(:, k), :); % Permuta envY de acordo com o resultado
        end
    end

    % Atualizando as matrizes, se global_corr for o último método
    if global_corr_last
        % Como é a última, não precisa atualizar muita coisa
        for k = find(~ind_conf)
            ind_conf(k) = true;
            P(:, k) = Ptmp(:, k);
            W(k, :, :) = W(k, P(:, k), :);
        end
        
        ind_solve{c_ind_solve} = ind_conf;
        c_ind_solve = c_ind_solve + 1;
    end
    
end

%% Correlação local
if must_do_localcorr
    if debug
        disp('--- Correlação local:')
    end

    converged = false;
    tmp = zeros(N, N); % Variável temporária para evitar o uso de squeeze
    perm = Ptmp;
    
    while(~converged)
 
        for k = find(~ind_conf)

            tmpRf = zeros(1, n_perms);
            for i = 1:n_perms

                tmptmp = zeros(N, N);
                
                % As frequências adjacentes são fixadas segundo a
                % permutação atual, e a frequência k em questão é permutada
                for kk =  1:length(Ycorr(k).indFreqs)
                    tmp(:, :) = Ycorr(k).AdjFreqCorr(kk, poss_perms(i, :), :);
                    tmptmp = tmptmp + tmp( :, perm(:, Ycorr(k).indFreqs(kk)) );
                end
                tmpRf(i) = sum(diag(tmptmp));

            end
            
            % Permutação que obteve maior correlação. Perceba que a
            % alteração na matriz perm de uma frequência também altera a
            % correlação das adjacentes com ela
            [dummy tmpind] = max(tmpRf);
            perm(:, k) = poss_perms(tmpind, :).';
            
        end

        % Testando convergência
        remaining_freqs = sum(sum(abs(perm(:, ~ind_conf) - Ptmp(:, ~ind_conf)), 1) ~= 0);
        if ~remaining_freqs
            converged = true;
        end

        % Atualizando as matrizes
        for k = find(~ind_conf)
            Ptmp(:, k) = perm(:, k); % Atualiza Ptmp
        end

    end
    
    % Atualizando as matrizes.
    for k = find(~ind_conf)
        ind_conf(k) = true;
        P(:, k) = Ptmp(:, k);
        W(k, :, :) = W(k, P(:, k), :);
        envY(k, :, :) = envY(k, P(:, k), :);
        Ymean(:, k) = Ymean(P(:, k), k);
        Yvar(:, k) = Yvar(P(:, k), k);
        Ycorr(k).AdjFreqCorr = Ycorr(k).AdjFreqCorr(:, P(:, k), :);
        
        % Atualiza a correlação das frequências adjacentes
        ind_changed = false(1, unique_K);
        ind_changed( getInvCustomInd(k, corr_neighbor_freq, 1, unique_K) ) = true;
        for kk = find(ind_changed)
            ind_adjchanged = find(Ycorr(kk).indFreqs == k);
            if isempty(ind_adjchanged),     continue;   end
            Ycorr(kk).AdjFreqCorr(ind_adjchanged, :, :) = Ycorr(kk).AdjFreqCorr(ind_adjchanged, :, P(:, k));
        end
    end
        
    ind_solve{c_ind_solve} = ind_conf;
    
end

%% Correlação com Threshold
if must_do_precorr
    if debug
        disp('--- Método Correlação de Frequências Adjacentes com Threshold:')
    end
    
    tmp = zeros(N, N); % Variável temporária para evitar o uso de squeeze
    ind_changed = true(1, unique_K); % Inicializa variável como se todas as permutações de todas as frequências tivessem mudado
    Rf = zeros(1, unique_K);
    
    % Enquanto houverem índices não confiáveis
    while ~isempty(find(~ind_conf, 1))

        % Encontra a correlação de todas as frequências de permutações não 
        % confiáveis com suas adjacentes. O lógico ind_changed otimiza
        % bastante o loop só encontrando a correlação em frequências que
        % foram modificadas e suas adjacentes
        for k = find(~ind_conf & ind_changed)
            
            tmpRf = zeros(1, n_perms);
            ind_tmpconf = ind_conf(Ycorr(k).indFreqs);

            for i = 1:n_perms
                tmp(:, :) = sum(Ycorr(k).AdjFreqCorr(ind_tmpconf, poss_perms(i, :), :), 1);
                tmpRf(i) = sum(diag(tmp));
            end
            
            [Rf(k) tmpind] = max(tmpRf);
            Ptmp(:, k) = poss_perms(tmpind, :).';
        end
        
        [dummy tmpind] = max(Rf);
        
        if dummy > corr_thre
            ind_conf(tmpind) = true;
            P(:, tmpind) = Ptmp(:, tmpind);
            Rf(tmpind) = 0;
            
            % Serve para otimização
            ind_changed = false(1, unique_K);
            ind_changed( getAdjInd(tmpind, corr_neighbor_freq, 1, unique_K) ) = true;

            % Ymean e Yvar não precisariam ser atualizadas, mas podem ser úteis
            % em alguma implementação futura
            W(tmpind, :, :) = W(tmpind, P(:, tmpind), :);
            envY(tmpind, :, :) = envY(tmpind, P(:, tmpind), :);
            Ymean(:, tmpind) = Ymean(P(:, tmpind), tmpind);
            Yvar(:, tmpind) = Yvar(P(:, tmpind), tmpind);
            Ycorr(tmpind).AdjFreqCorr = Ycorr(tmpind).AdjFreqCorr(:, P(:, tmpind), :);
            
            % Atualiza a correlação das frequências adjacentes
            for kk = find(ind_changed)
                ind_adjchanged = find(Ycorr(kk).indFreqs == tmpind);
                if isempty(ind_adjchanged),     continue;   end
                Ycorr(kk).AdjFreqCorr(ind_adjchanged, :, :) = Ycorr(kk).AdjFreqCorr(ind_adjchanged, :, P(:, tmpind));
            end
        else
            break;
        end
    end

    ind_solve{c_ind_solve} = ind_conf;
    c_ind_solve = c_ind_solve + 1;
    
    if debug
        disp( sprintf('Número de bins resolvidos: %.0f\n', sum(ind_conf)) )
    end
end

%% Correlação Harmônica
if must_do_harmonic
    if debug
        disp('--- Método Correlação de Frequências Harmônicas:')
    end

    env_f = zeros(N, num_of_frames);
    
    % Encontra a correlação de todas as frequências de permutações não confiáveis com
    % suas harmônicas
    for k = find(~ind_conf)
        env_f(:, :) = envY(k, :, :);
        
        ind_harm = false(size(ind_conf));
        if(k <= K/2 + 1) % Se a frequência for positiva ou -fs/2
            ind_harm(getCustomHarmInd(k, 1, K/2 + 1)) = true;
        else
            ind_harm(getCustomHarmInd(k, K + 1, K/2 + 1)) = true;
        end
        env_g = envY(ind_harm & ind_conf, :, :);

        % Corre as frequências adjacentes encontrando as correlações
        if isempty(env_g)
            continue; % Se não houver harmônicos com permutação confiável
        else
            [Rf Ptmp(:, k)] = maxcorr(Ytmp, env_g);
        end

        if Rf > harmonic_thre
            ind_conf(k) = true;
            P(:, k) = Ptmp(:, k);
            W(k, :, :) = W(k, P(:, k), :);
            envY(k, :, :) = envY(k, P(:, k), :);
            Ymean(:, k) = Ymean(P(:, k), k);
            Yvar(:, k) = Yvar(P(:, k), k);
            Ycorr(k).AdjFreqCorr = Ycorr(k).AdjFreqCorr(:, P(:, k), :);
            
            % Atualiza a correlação das frequências adjacentes            
            ind_changed = false(1, unique_K);
            ind_changed( getAdjInd(tmpind, corr_neighbor_freq, 1, unique_K) ) = true;
            for kk = find(ind_changed)
                ind_adjchanged = find(Ycorr(kk).indFreqs == k);
                if isempty(ind_adjchanged),     continue;   end
                Ycorr(kk).AdjFreqCorr(ind_adjchanged, :, :) = Ycorr(kk).AdjFreqCorr(ind_adjchanged, :, P(:, k));
            end
        
        end
    end

    ind_solve{c_ind_solve} = ind_conf;
    c_ind_solve = c_ind_solve + 1;

    if debug
        disp( sprintf('Número de bins resolvidos: %.0f\n', sum(ind_conf)) )
    end
end

%% Correlação sem Threshold
if must_do_corr 
    if debug
        disp('--- Método Correlação de Frequências Adjacentes sem Threshold:')
    end
    
    tmp = zeros(N, N); % Variável temporária para evitar o uso de squeeze
    ind_changed = true(1, unique_K); % Inicializa variável como se todas as permutações de todas as frequências tivessem mudado
    Rf = zeros(1, unique_K);
    
    % Enquanto houverem índices não confiáveis
    while ~isempty(find(~ind_conf, 1))

        % Encontra a correlação de todas as frequências de permutações não 
        % confiáveis com suas adjacentes. O lógico ind_changed otimiza
        % bastante o loop só encontrando a correlação em frequências que
        % foram modificadas e suas adjacentes
        for k = find(~ind_conf & ind_changed)
            
            tmpRf = zeros(1, n_perms);
            ind_tmpconf = ind_conf(Ycorr(k).indFreqs);

            if (corr_conf_condition) % Se somente devem ser calculadas as somas de correlações para frequências com permutação já resolvida...
                for i = 1:n_perms
                    tmp(:, :) = sum(Ycorr(k).AdjFreqCorr(ind_tmpconf, poss_perms(i, :), :), 1);
                    tmpRf(i) = sum(diag(tmp));
                end
            else % ... ou não
                for i = 1:n_perms
                    tmp(:, :) = sum(Ycorr(k).AdjFreqCorr(:, poss_perms(i, :), :), 1);
                    tmpRf(i) = sum(diag(tmp));
                end
            end
            
            [Rf(k) tmpind] = max(tmpRf);
            Ptmp(:, k) = poss_perms(tmpind, :).';
        end
        
        [dummy tmpind] = max(Rf);
        ind_conf(tmpind) = true;
        P(:, tmpind) = Ptmp(:, tmpind);
        Rf(tmpind) = 0;
        
        % Serve para otimização
        ind_changed = false(1, unique_K);
        ind_changed( getInvCustomInd(tmpind, corr_neighbor_freq, 1, unique_K) ) = true;
        
        % Ymean e Yvar não precisariam ser atualizadas, mas podem ser úteis
        % em alguma implementação futura
        W(tmpind, :, :) = W(tmpind, P(:, tmpind), :);
        envY(tmpind, :, :) = envY(tmpind, P(:, tmpind), :);
        Ymean(:, tmpind) = Ymean(P(:, tmpind), tmpind);
        Yvar(:, tmpind) = Yvar(P(:, tmpind), tmpind);
        Ycorr(tmpind).AdjFreqCorr = Ycorr(tmpind).AdjFreqCorr(:, P(:, tmpind), :);
        
        % Atualiza a correlação das frequências adjacentes
        for kk = find(ind_changed)
            ind_adjchanged = find(Ycorr(kk).indFreqs == tmpind);
            if isempty(ind_adjchanged),     continue;   end
            Ycorr(kk).AdjFreqCorr(ind_adjchanged, :, :) = Ycorr(kk).AdjFreqCorr(ind_adjchanged, :, P(:, tmpind));
        end

    end
    
    ind_solve{c_ind_solve} = ind_conf;

end

%% Supervisionado
if must_do_supervised
    if debug
        disp('--- Método Supervisionado (by Victorio):')
    end
    
    for k = find(~ind_conf)

        % Compara o sinal encontrado com o original para tentar descobrir
        % se houve permutação ou não. A maior correlação corresponde à
        % permutação correta
        env_f = envelope( squeeze(Y(k, :, :)), 1 );
        env_g(1, :, :) = envelope( squeeze(S(k, :, :)), 1 );
        [dummy P(:, k)] = maxcorr(env_f, env_g);
        
    end
end

%% MAXSIR (do Makino)
if must_do_maxSIR
    if debug
        disp('--- Método MaxSIR (Supervisionado by Makino):')
    end
    
    for k = find(~ind_conf)

        % Utiliza o método de Makino de encontrar o traço máximo de uma
        % função (perm*W*Q).^2
        [dummy P(:, k)] = maxSIR(squeeze(W(k, :, :)), squeeze(Q(k, :, :)));
        
    end
end

%% Cálculos finais
% Por causa da simetria da FFT (se for utilizada)
for k = unique_K+1:K
    P(:,k) = P(:,K+2-k);
end

% Organizando a saída ind_solve
for i = length(ind_solve):-1:2
    ind_solve{i} = xor(ind_solve{i}, ind_solve{i-1});
end
    
clear global SPEED_OF_SOUND

%% Functions obsoletas

%{
function indArray = getHarmInd(f, N)
% indArray = getHarmInd(f, N) retorna os índices dos harmônicos
% da frequência cujo índice é 'f'. Considera-se que o primeiro elemento do
% vetor de frequências corresponde à frequência 0, e o segundo, a fs/'N',
% onde fs é a frequência de amostragem e 'N' o número de frequências (ou
% de bins da FFT), o terceiro, a (2*fs)/'N' e assim por diante.
%
% A saída é um vetor LINHA
%
% Exemplos:
% getHarmInd(1, 16) retorna o vetor [], pois o índice 1 corresponde à
%                   frequência 0, que não possui harmônicos
% getHarmInd(2, 16) retorna o vetor [3 4 5 6 7 8 9 10 11 12 13 14 15 16], que
%                  corresponde ao vetor inteiro
% getHarmInd(5, 16) retorna o vetor [9 13]


if f == 1
    indArray = [];
else
     % O número de harmônicos que uma dada frequência possui é o número de bins entre a frequência máxima e o bin específico (N-f),
     % dividido pelo tamanho do pulo do harmônico, que no caso, é igual ao índice da frequência (f-1, pois quando f é 1 o índice é 0)
    num_harm = floor( (N-f) / (f-1) );
    harmArray = 2:num_harm+1;
    indArray = (f - 1) * harmArray + ones(1, length(harmArray)); % Os harmônicos são (f-1)*2 + 1, (f-1)*3 + 1, (f-1)*4 + 1, e assim por diante.
end
%}