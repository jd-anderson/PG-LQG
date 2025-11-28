function grad = gradient_modelfree(A, B, C, Q, R, K, p, r, n_s, V, W)
%GRADIENT_MODELFREE Zero-order gradient estimate aligned with text
%   grad = gradient_modelfree(A,B,C,Q,R,K,p,r,n_s,V,W)
% Inputs:
%   A,B,C, Q,R     - system and cost matrices
%   K              - current lifted gain (n_u × (p*(n_u+n_y)))
%   p              - history length
%   r              - smoothing radius (Frobenius norm)
%   n_s            - number of samples for ZO estimator
%   V, W           - measurement/process noise covariances for cost
%
% Output:
%   grad           - estimated gradient, same size as K

% Determine sizes
[m, n] = size(K);     % m = n_u, n = p*(n_u+n_y)
n_u = m;
n_y = size(C,1);

dif = zeros(m, n);

% Optional: baseline subtraction for variance reduction
% J0 = cost_modelfree(A,B,C,K,Q,R,p,V,W);

for i = 1:n_s
    % sample U on Frobenius sphere radius r
    U = randn(m, n);
    U = U / norm(U, 'fro') * r;
    
    % evaluate cost at K+U
    Jp = cost_modelfree(A, B, C, K + U, Q, R, p, V, W);
    % accumulate one-point estimator
    dif = dif + Jp * U;
    % with baseline: dif = dif + (Jp - J0) * U;
end

% ZO factor: p * n_u * (n_u + n_y) / (n_s * r^2)
factor = p * n_u * (n_u + n_y) / (n_s * r^2);
grad = factor * dif;

% Clip gradient if too large (you can choose threshold)
max_grad = 1e3;  % for example; adjust as needed
ng = norm(grad, 'fro');
if ng > max_grad
    grad = grad * (max_grad / ng);
    a =1;
end

end
