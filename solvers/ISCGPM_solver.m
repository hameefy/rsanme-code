function [x, info] = ISCGPM_solver(u0, X, Y, eta, params)
maxIter = params.maxIter; tol = params.tol; sigma = params.sigma; rho = params.rho;
alphaBar = params.alphaBar; nu = params.nu; mu = params.mu; r = params.r; gamma = params.gamma; etaSpec = params.etaSpec; %#ok

x_prev = u0; x = u0; k = 0; resHist = zeros(maxIter+1,1);
fevals = 0;                                    % initialize function-eval counter

F_y_prev = logistic_residual(X, Y, x_prev, eta); fevals = fevals + 1;
resHist(1) = norm(F_y_prev); d_prev = -F_y_prev;
eta_sp = 0.3; nu_sp = nu;

while k < maxIter
    if k == 0 || norm(x - x_prev) == 0
        alpha_k = alphaBar;
    else
        alpha_k = min(alphaBar, norm(x - x_prev)^2 / (k^2));
    end
    y = x + alpha_k * (x - x_prev);
    F_y = logistic_residual(X, Y, y, eta); fevals = fevals + 1;
    if norm(F_y) <= tol, x = y; break; end

    if k == 0
        d = -F_y; theta_k = 1; beta_k = 0; %#ok
    else
        s = y - (x_prev + alpha_k * (x_prev - x_prev));
        w = F_y - F_y_prev + r * s;
        theta_hat = (s' * w) / max(s' * s, 1e-16);
        theta_k = min(max(theta_hat, eta_sp), nu_sp);
        beta_k = (F_y' * w) / max(max(mu * norm(d_prev) * norm(w), norm(F_y)^2), 1e-16);
        d = -theta_k * F_y + beta_k * d_prev;
    end
    t = 1;
    while true
        z = y + t * d;
        F_z = logistic_residual(X, Y, z, eta); fevals = fevals + 1;
        lhs = -F_z' * d; rhs = sigma * t * norm(F_z) * (norm(d)^2);
        if lhs >= rhs, break; else, t = rho * t; end
        if t < 1e-16, break; end
    end
    xi = (F_z' * (y - z)) / max(norm(F_z)^2, 1e-16);
    x_next = y - gamma * xi * F_z;
    k = k + 1; resHist(k+1) = norm(F_z);
    x_prev = x; x = x_next; F_y_prev = F_y; d_prev = d;
    if norm(F_z) <= tol, break; end
end

info.iter = k;
info.fval = fevals;            % number of residual evaluations
info.resHist = resHist(1:k+1); % history of the norms
info.resk = norm(F_z); % the final norm

