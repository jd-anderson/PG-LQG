function S = calculate_S(A, B, C, Q, R, W, V, p)
% CALCULATE_S  Compute the history-to-state mapping S* for LQG (history representation)
%
%   S = calculate_S(A, B, C, Q, R, W, V, p)
%
% Inputs:
%   A, B, C   - system matrices (A: nx×nx, B: nx×nu, C: ny×nx)
%   Q, R      - cost weights (for LQR in estimator dynamics)
%   W, V      - process and measurement noise covariances (nx×nx, ny×ny)
%   p         - history length (integer, ≥ observability index)
%
% Output:
%   S         - mapping matrix S* of size nx × [p*(nu + ny)], such that
%               x̂_t = S * z_{t,p}, where z_{t,p} = [u_{t-1};…;u_{t-p}; y_{t-1};…;y_{t-p}].
%
% This implements the formula in the history representation (see :contentReference[oaicite:0]{index=0}):
%   Ã = (I - L C) A,   B̃ = (I - L C) B,   where L is the steady-state Kalman gain.
%   Compute L via discrete-time Riccati: Σ = DARE(A', C', W, V); L = Σ C' (C Σ C' + V)^{-1}.
%   Compute estimator LQR gain K_star via dlqr on (Ã, B̃) with cost Q̃ = C'Q C, R.
%   Build blocks:
%     Fu,p = [B̃, Ã B̃, …, Ã^(p-1) B̃]   (nx × p*nu)
%     Fy,p = [L, Ã L, …, Ã^(p-1) L]     (nx × p*ny)
%     Tu,p, Ty,p, Ox,p per paper (see below).
%   Then S = [Fu,p + Ã^p * Ox,p^† * (I - Tu,p),   Fy,p - Ã^p * Ox,p^† * Ty,p].
%
% Example:
%   S = calculate_S(A, B, C, Q, R, W, V, 3);
%
% Reference: On the Gradient Domination of the LQG Problem, Sec. II.A :contentReference[oaicite:1]{index=1}.

    %% Dimensions
    nx = size(A,1);
    nu = size(B,2);
    ny = size(C,1);

    %% 1. Compute steady-state Kalman gain L
    % Solve discrete-time Riccati for estimator: A', C', W, V
    Sigma_tilde = dare(A', C', W, V);
    L = Sigma_tilde * C' / (C * Sigma_tilde * C' + V);  % nx×ny

    %% 2. Estimator dynamics
    I_nx = eye(nx);
    A_tilde = (I_nx - L*C) * A;  % nx×nx
    B_tilde = (I_nx - L*C) * B;  % nx×nu

    %% 3. Compute LQR gain K_star for estimator dynamics
    % Cost on estimated state: Q̃ = C'Q C
    Q_tilde = C' * Q * C;
    [K_star, ~, ~] = dlqr(A_tilde, B_tilde, Q_tilde, R);  % returns K_star so u = -K_star * x̂

    %% 4. Build Fu,p and Fy,p
    Fu = zeros(nx, p*nu);
    Fy = zeros(nx, p*ny);
    A_power = eye(nx);
    for i = 1:p
        % Fu(:, (i-1)*nu+1 : i*nu) = A_tilde^(i-1) * B_tilde
        Fu(:, (i-1)*nu+1 : i*nu) = A_power * B_tilde;
        % Fy(:, (i-1)*ny+1 : i*ny) = A_tilde^(i-1) * L
        Fy(:, (i-1)*ny+1 : i*ny) = A_power * L;
        % update A_power = A_tilde^i for next iteration
        A_power = A_tilde * A_power;
    end

    %% 5. Build Tu,p and Ty,p
    % Tu,p is (p*nu × p*nu), Ty,p is (p*nu × p*ny)
    Tu = zeros(p*nu);
    Ty = zeros(p*nu, p*ny);
    % According to paper: entries at block (row, col) when col == row+1:
    %   Tu_block = K_star * B_tilde, Ty_block = K_star * L
    for row = 1:p
        col = row + 1;
        if col <= p
            row_idx = (row-1)*nu + (1:nu);
            col_idx_u = (col-1)*nu + (1:nu);
            col_idx_y = (col-1)*ny + (1:ny);
            Tu(row_idx, col_idx_u) = K_star * B_tilde;  % nu×nu
            Ty(row_idx, col_idx_y) = K_star * L;        % nu×ny
        end
    end

    %% 6. Build Ox,p
    % Ox is (p*nu × nx). According to corrected assignment:
    %   For i = 1..p: block rows (i-1)*nu+1 : i*nu equals K_star * A_tilde^(p-i)
    Ox = zeros(p*nu, nx);
    for i = 1:p
        % power = A_tilde^(p-i)
        A_pow = A_tilde^(p - i);
        Ox((i-1)*nu+1 : i*nu, :) = K_star * A_pow;  % (nu×nx)
    end

    %% 7. Compute pseudoinverse Ox†
    % Use pinv for numerical robustness: Ox_pinv is nx×(p*nu)
    Ox_pinv = pinv(Ox);

    %% 8. Compute A_tilde^p
    A_tilde_p = A_tilde^p;

    %% 9. Compute S* blocks
    I_pnu = eye(p*nu);
    % M1 = Fu + Ã^p * Ox† * (I - Tu)
    M1 = Fu + A_tilde_p * Ox_pinv * (I_pnu - Tu);
    % M2 = Fy - Ã^p * Ox† * Ty
    M2 = Fy - A_tilde_p * Ox_pinv * Ty;

    %% 10. Form S*
    S = [M1, M2];  % size nx × [p*nu + p*ny] = nx × p*(nu+ny)
end

