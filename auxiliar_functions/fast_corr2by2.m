function R = fast_corr2by2(P, Q, mP, mQ, vP, vQ)
%
% R = fast_corr2by2(P, Q, mP, mQ, vP, vQ) encontra correlações entre vetores
% das matrizes P e Q
%
% P, Q - matrizes com um vetor por linha, ou seja,
%   P = [p1; p2; p3 ... pm], q = [q1; q2; q3 ... qn]
% mP, mQ, vP, vQ - vetores coluna com as médias e variâncias de cada vetor
% 
% R - matriz m x n de correlação da forma: 
%   R = [corr(p1,q1) corr(p1,q2) corr(p1,q3) ... corr(p1,qn)
%        corr(p2,q1) corr(p2,q2) corr(p2,q3) ... corr(p2,qn)
%                              ...
%        corr(pm,q1) corr(pm,q2) corr(pm,q3) ... corr(pm,qn)]
%
m = size(P, 1);
n = size(Q, 1);
N = size(P, 2);

P = conj(P);
mP = conj(mP);
R = zeros(m, n);
for cn = 1:n
    R(:, cn) = (sum(P .* repmat(Q(cn, :), m, 1), 2)/(N-1) - (N/(N-1)) * mP * mQ(cn)) ./ (sqrt(vP * vQ(cn)));
end