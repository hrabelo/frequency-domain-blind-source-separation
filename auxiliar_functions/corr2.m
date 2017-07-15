function r = corr2(x, y)
%
% r = corr2(x, y) encontra o coeficiente de correlação entre x e y
%

x = x(:).';
y = y(:).';
mx= mean(x);
my = mean(y);
vx = var(x);
vy = var(y);

r = fast_corr2by2(x, y, mx, my, vx, vy);