function [z,V] = pre_whitening(x)
%
% Função que efetua um branqueamento dos dados via decomposição por autovalores
% Utilize a função pre_centering antes dessa, senão a média dos valores não
% será nula.
%
% Sintaxe:
%  [z,V] = pre_whitening(x)
% 
% Argumento de entrada:
%  x -> dados (possivelmente complexos, dispostos a cada linha) a serem branqueados 
% 
% Argumento de saída:
%  z -> dados branqueados
%  V -> matriz que efetua a transformação linear que branqueia x (z = V * x)

% Última modificação: 24/01/10 - corrigido bug do transpose
%                     07/07/10 - agora a média do sinal branqueado é 0,
%                                utilizando a função pre_centering

%% Inicialização
N = size(x,2);
M = size(x,1);
R = zeros(M);

%% Cálculo da matriz de covariância
for ind = 1:N
    R = R + x(:,ind) * x(:,ind)'; % O operador ' já representa transposição hermitiana, ou seja, usar também a função conj é errado
end
 
R = R/N;

%% Autovalores (D) e Autovetores (E)
[E,D] = eig(R);
 
% Pondo em ordem decrescente os autovalores (ATENÇÃO: supõe-se que não haja autovalores idênticos)

[d ind] = sort(diag(D), 'descend');
D = diag(d);
E = E(:,ind);


%% Calculando a matriz de transformação
V = D^(-1/2) * E';
z = V * x;
