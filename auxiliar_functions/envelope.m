function env = envelope(y, type, w)
%
% env = envelope(y, type) calcula o envelope de y, que é um ou mais vetores,
% dispostos a cada linha. NÃO FUNCIONA COM MAIS DE 2 dimensões. O parâmetro
% 'type' diz o tipo de envelope. 1 é valor absoluto e 2 é valor relativo
% quadrático em relação aos outros elementos da mesma coluna.
%
% env = envelope(y, type, w), é para o método 3, onde w é uma matriz com o
% mesmo número de linhas que y, e um número de colunas igual ao número de
% sensores.
%
%   1 ('AbsValue')  - Valor absoluto.
%   2 ('PowValue')  - Razão entre a potência do valor do vetor e a soma das
%                     potências do mesmo valor em todos os vetores. Por
%                     definição, varia entre 0 e 1.
%   3 ('PowValue2') - Similar ao acima, porém cada elemento j de uma linha
%                     i de y, onde a = w^(-1) e M é o número de sensores, é
%                     sum(||a(1,i)*y(i,j) + a(2,i)*y(i,j) + ... + a(M,i)*y(i,j)||^2)
%   4 ('SPDValue')   - 

% É aconselhável que esta função não dependa da fase do sinal.

N = size(y, 1); % Número de vetores
env = zeros(size(y));

if type == 1
    env = abs(y);

elseif type == 2
    for i = 1:N
        env(i,:) = y(i, :).*conj(y(i, :));
    end
    env = bsxfun(@rdivide, env, sum(env, 1));

    % Pre MatLab 2008
    %sumenv = repmat(sum(env, 1), size(env, 1), 1);
    %env = env ./ sumenv;
    
elseif type == 3
    a = inv(w);
    for i = 1:N
        tmpy = a(:,i) * y(i,:);
        env(i, :) = sum(tmpy.*conj(tmpy), 1);
    end
    %env = bsxfun(@rdivide, env, sum(env, 1));
    
    % Pre MatLab 2008
    sumenv = repmat(sum(env, 1), size(env, 1), 1);
    env = env ./ sumenv;
    
elseif type == 4
    for i = 1:N
        env(i, :) = y(i, :).*conj(y(i, :));
    end
    env = 10*log(env);
%    env = bsxfun(@rdivide, env, sum(env, 1));
    
    % Pre MatLab 2008
    %sumenv = repmat(sum(env, 1), size(env, 1), 1);
    %env = env ./ sumenv;

end