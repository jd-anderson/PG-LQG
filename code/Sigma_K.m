function Sigma_nu = Sigma_K(A, C, W, V)
% COMPUTE_SIGMA_NU  Compute Sigma_nu = L*(C*Sigma_tilde*C' + V)*L'
% for the Kalman estimator, given system matrices.
%
%   Sigma_nu = compute_Sigma_nu(A, C, W, V)
%
% Inputs:
%   A   - state transition matrix (nx×nx)
%   C   - output matrix (ny×nx)
%   W   - process noise covariance (nx×nx)
%   V   - measurement noise covariance (ny×ny)
%
% Output:
%   Sigma_nu - covariance of innovation term: L*(C*Sigma_tilde*C' + V)*L'
%
% Steps:
%   1. Solve estimator Riccati: Sigma_tilde = dare(A', C', W, V)
%   2. Compute Kalman gain: L = Sigma_tilde*C'/(C*Sigma_tilde*C' + V)
%   3. Sigma_nu = L*(C*Sigma_tilde*C' + V)*L'
%
% Example:
%   Sigma_nu = compute_Sigma_nu(A, C, 0.01*eye(size(A)), 0.01*eye(size(C,1)));
%
    % Check dimensions
    nx = size(A,1);
    ny = size(C,1);
    assert(all(size(A) == [nx, nx]), 'A must be square nx×nx');
    assert(size(C,2) == nx, 'C must be ny×nx');
    assert(all(size(W) == [nx, nx]), 'W must be nx×nx');
    assert(all(size(V) == [ny, ny]), 'V must be ny×ny');

    % 1. Solve estimator Riccati: Sigma_tilde
    Sigma_tilde = dare(A', C', W, V);  % solves A' X A - X - A' X C'*(C X C' + V)^{-1} C X A + W = 0

    % 2. Kalman gain
    S = C * Sigma_tilde * C' + V;  % ny×ny
    L = Sigma_tilde * C' / S;      % nx×ny

    % 3. Sigma_nu
    Sigma_nu = L * S * L';         % nx×nx
end
