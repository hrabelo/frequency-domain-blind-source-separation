function checkwin(wind, hop, span)
% checkwin(wind, hop) plota um gráfico para checar se a janela "wind" atende
% às especificações de Constant-Overlap-Add quando utilizada com um pulo
% "hop"
%
% Última modificação: 24/08/2010

M = length(wind);
if nargin < 3
    span = 5*M;                 % overlap-add span
end

z = zeros(span,1);  figure;  hold on;  s = z;
for so=0:hop:span-M
  ndx = so+1:so+M;              % current window location
  s(ndx) = s(ndx) + wind;       % window overlap-add
  wzp = z; wzp(ndx) = wind;     % for plot only 
  plot(wzp,'--k');              % plot just this window
end
plot(s,'k');  hold off;  % plot window overlap-add
%legend('A linha azul representa o Overlap-Add', 'FontSize', 12)
xlim([1 span]); xlabel('Amostras');
ylim([0 1.1*max(s)]);