function [y, u, evol_u, nit] = cunICA(z, type, Max_It)
% COMPLEX UNITARY ICA:  [y, u, evol_u] = cunICA(x)
%                       [y, u, evol_u] = cunICA(x, type, Max_It)
% u -> matriz de separação unitária
% y -> fontes separadas
% evol_u -> evolução da matriz de separação com o tempo, as iterações estão
% na 3ª dimensão
% nit -> número de iterações para convergir
%
% z -> matriz de misturas (cada linha é uma mistura)
% type -> função G a utilizar. 1 - genLaplace
% Max_It -> numero máximo de iterações
%
% SOMENTE RESOLVE , POR ENQUANTO, COM O MESMO NÚMERO DE FONTES QUE MISTURAS

%% Processando as entradas e saídas
if nargout > 2,     save_progress = 1;
else                save_progress = 0;                              end
    
if nargin < 3
    Max_It = 250;
    if nargin < 2
        type = 1;
    end
end

%% Inicialização
N = size(z, 1);
num_of_samples = size(z, 2);

u = eye(N);
y = z;
alpha = 0.1;

convergence_test = true;
Min_It = 8; % Número de iterações a partir do qual se testa a convergência
err_thre = 0;  % Threshold do erro para teste de convergência (default 0)
err_u = zeros(N,1);

% A densidade estimada das fontes é p(x) = C*exp( -(sqrt(abs(x)^2 + alpha) / b) )
% C não é importante, e b é a variância, que somente afeta a escala do
% sinal, que é ajustada depois em um FDBSS.
% Portanto, a densidade fica
%       p(x) = exp(-sqrt(abs(x)^2 + alpha))

if save_progress,       evol_u = zeros(Max_It, N, N);               end  % LIMITA NUMERO DE FONTES IGUAL NUMERO DE MISTURAS

%% Fast ICA
for it = 1:Max_It
    mod_y2 = y.*conj(y);
    ant_u = u;
    
    switch type
        case 1
        % Encontrando as funções g(y) e g'(y)
        tmp = 0.5 ./ sqrt(mod_y2 + alpha);
        phi = tmp .* y; % na verdade é tmp .* conj(y), mas fazendo desta forma não precisamos encontrar phi(y)
        gg = tmp .* (1 - 0.5*mod_y2 ./ (mod_y2 + alpha));

        % Encontrando a matriz gamma e phi(y), que é igual conj(g), daí se
        % torna desncessário
        gamma = diag( sum(gg, 2) / num_of_samples );
    end
    
    % Atualizando a matriz de separação u
    u = gamma * u - phi * z' / num_of_samples;
    
    % Ajustando para uma matriz unitária
    u = (u*u')^(-1/2) * u;
    
    if convergence_test
        %% Descobrindo a convergência
        % Basta testar se o produto interno do vetor atual e o antigo é 1,
        % pois isto significa que eles apontam na mesma direção
        for cn = 1:N
            err_u(cn) = abs(u(cn, :) * ant_u(cn, :)');
        end
        
        % Como cada vetor é unitário, o produto interno não pode ser maior
        % que 1. Se não existir mais nenhum valor menor do que 1, sai do
        % loop. Espero um número mínimo de iterações.
        if it > Min_It
            if ~sum(err_u < (1 - err_thre)),                  break;      end
        end
        
    end
    
    if save_progress,       evol_u(it, :, :) = u;                   end
    
    y = u * z;
end

nit = it;
