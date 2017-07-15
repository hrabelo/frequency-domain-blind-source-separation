%% PATH
Path = 'C:\Users\Hian\Documents\final-project\Matlab\commons\audio_files\';

% Cria arrays 1x2;
Sources = cell(1, 2);

% Obtém os arquivos de áudi
Sources{1} = [Path 'female_src_1.wav'];
Sources{2} = [Path 'male_src_1.wav'];

%% PARÂMETROS DA STFT

K = 4096;                                   % Número de bins da FFT utilizada (pontos)
windlen = 2048;                             % Comprimento da janela

Jc = cell(1,5);
Jc{1} = windlen / 2;
Jc{2} = windlen / 4;                        % Pulo da STFT
Jc{3} = windlen / 4;
Jc{4} = windlen / 4;
Jc{5} = windlen / 4;
Jc{6} = windlen;
Jc{7} = windlen / 2;
Jc{8} = windlen / 4;
J = Jc{8};

wind = cell(1,5);
wind{1} = hann(windlen, 'periodic').';      % Janelas que atendem à COLA se J = len/2, onde len é o comprimento da janela.
wind{2} = 0.5*hann(windlen, 'periodic').';  % Janelas que atendem quando J = len/4
wind{3} = 0.677*chebwin(windlen).';
wind{4} = 0.72*window(@blackmanharris, windlen).';
wind{5} = 0.71*window(@nuttallwin, windlen).';
wind{6} = ones(1,windlen);
wind{7} = 0.5*ones(1,windlen);
wind{8} = 0.25*ones(1,windlen);

winda = wind{8};                            % Janela de análise
winds = ones(1,windlen);                    % Janela de síntese
wdft_par = 0;                               % Se utilizada WDFT, coloque um valor diferente de 0, mas entre -1 e 1 exclusive.
debug = 1;                                  % Depuração da rotina BSS_FD

%% PARÂMETROS DO ICA
sep_meth = 'natica';                         % Método ICA utilizado
Num_It = 250;                               % Número máximo de iterações do ICA
eta = .2;                                   % Passo de adaptação do ICA
natica_meth = 'sign';                       % Parâmetro do natica e conjica
nonholonomic = false;                       % Parâmetro do natica e conjica
laplace_alpha = 0.1;                        % Parâmetro do natica_laplace
sourcepdf_dev = 1;                          % Parâmetro do natica (desvio padrão aproximado da pdf do sinal da fonte)


%% PARÂMETROS PARA A MATRIZ DE MISTURAS, ENTRADA PARA O SCRIPT BSS_READ
mixing_mode = 'ismbuild';                   % Método para achar as misturas. Ver script bss_read para informações
N = 2;                                      % Número de fontes
M = 2;                                      % Número de misturas
num_samp = 160000;                           % Número de amostras a ler
reverb_time = 0.1;                         % Tempo de reverberação
mics_struct = 'cluster2d';                  % Estrutura de montagem dos microfones
dist_scrmic = 1;                            % Distância de cada fonte até os microfones, em metros
ang_src(1) = 45;                           % Se utilizado ISM, ângulo da fonte 1, em graus
ang_src(2) = 120;                           % Se utilizado ISM, ângulo da fonte 2, em graus
Source{1} = Sources{1};                     % Arquivo da fonte 1 
Source{2} = Sources{2};                     % Arquivo da fonte 2
plotroom_flag = 1;

%% PARÂMETROS PARA RESOLVER A PERMUTAÇÃO
perm_meth = 'tdoa';               % Método para resolver o problema. Digite help bss_fdpermsolve 
dist_mics = 0.04;                           % Distância entre os microfones, em metros (também é ENTRADA para o BSS_READ, se se utilizar simulação)
                                            % Distância entre os microfones das pontas
%dist_max_mics = 0.04;                       % É saída da função BSS_READ quando mixing_mode = 'ismbuild'. Anote quando rodar o build
dirpat_thre = 0;                          % Uma das condições do método DOA. Melhor deixar -Inf e não utilizá-la
corr_thre = 0.3*N*6;                        % Threshold que diz se o método de correlação entre frequências adjacentes é confiável. O valor máximo é N*2*NumAdjFrequencies, com N=2, é 12. Usei 30% do máximo
                                            % NumAdjFrequencies é um parâmetro da função bss_fdpermsolve, que deixei como 'default, ou seja, 3
corr_env = 'PowValue2';                      % Envelope para calcular a correlação ('AbsValue' ou 'PowValue')
harm_thre = 0.1*N*6;                        % Threshold que diz se o método de correlação entre harmônicos é confiável. O valor máximo é N*6. Usei 10% do valor máximo

%% SUAVIZANDO OS FILTROS
smoothFlag = 1;                             % Se deve suavizar os filtros
%smooth_filter = [0.25 0.5 0.25]; % Hanning window
smooth_filter = [0.003 0.0602 0.2516 0.3902 0.2516 0.0602 0.003]; % Chebyshev window
%smooth_filter = [0.01 0.0817 0.24 0.3363 0.24 0.0817 0.01]; % Blackman window
%smooth_filter = [0.0092 0.0795 0.2407 0.3409 0.2407 0.0795 0.0092]; % Nuttall window
%smooth_filter = [0.0014 0.0032 0.0129 0.9787 0.0129 0.0032 0.0014]; % Kaiser window
                                            % Filtro utilizado (padrão [0.25 0.5 0.25])

%% EXECUTA A ROTINA

disp(perm_meth)

source_combs = nchoosek(1:8, N);

BSS_FD;

[SDR, SIRb, SAR, perm]=bss_eval_sources(x,s);
[SDR, SIR, SAR, perm]=bss_eval_sources(y,s);   disp(SIR - SIRb); disp(SDR); disp(SAR)

%[SDR,SIR,SAR,perm]=bss_eval_sources(yconv(:,1:48000),s);   disp(SIR); disp(SDR); disp(SAR)
%[SDR,SIR,SAR,perm]=bss_eval_sources(yconvt(:, 1:48000),s);   disp(SIR); disp(SDR); disp(SAR)

clear tmpc

%% RASCUNHO
%doa_acertou = metodo_acertou & ind_solve{1};        disp(sum(doa_acertou)/sum(ind_solve{1}))
%precorr_acertou = metodo_acertou & ind_solve{2};    disp(sum(precorr_acertou)/sum(ind_solve{2}))
%harm_acertou = metodo_acertou & ind_solve{3};       disp(sum(harm_acertou)/sum(ind_solve{3}))
%corr_acertou = metodo_acertou & ind_solve{4};       disp(sum(corr_acertou)/sum(ind_solve{4}))
