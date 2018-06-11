function [Gs_stat_mu, P] = save_gs_stat(Y)

d = size(Y, 2)-1;
N = size(Y, 1);

w = Y(:, d+1);
X = Y(:, 1:d);

Gs_stat_mu = w' * X;

Gs_stat_R = bsxfun(@minus, X, Gs_stat_mu);

% p_ind = residualR_m((1:N)', w, N);

p_ind = randsample((1:N)', N, true, w);
Gs_stat_R = Gs_stat_R(p_ind, :);

P = Gs_stat_R' * Gs_stat_R/ size(X, 1);

min_eig = min(eig(P));
if min_eig <= 1e-15
    P = P + (-min_eig + 1e-15) * eye(size(Gs_stat_R,2));
end