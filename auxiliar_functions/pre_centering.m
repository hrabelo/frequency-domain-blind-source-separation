function x = pre_centering(x)
%
% Centraliza os dados (torna sua média zero)
%
% Sintaxe:
%  x = pre_centering(x)
% 
% Argumento de entrada:
%  x -> dados (possivelmente complexos, dispostos a cada linha) a serem branqueados 
% 
% Argumento de saída:
%  x -> dados centralizados
%

%x = bsxfun(@minus, x, mean(x, 2));

% Pre MatLab 2008
medx = repmat(mean(x, 2), 1, size(x, 2));
x = x - medx;
