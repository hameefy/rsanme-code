function [x, info] = liu2022(F, x0, opts)
%LIU2022  Algorithm 2.1 of J.K. Liu, Z.L. Lu, J.L. Xu, S. Wu, Z. Tu, "An efficient
%         projection-based algorithm without Lipschitz continuity for large-scale
%         nonlinear pseudo-monotone equations", J. Comput. Appl. Math. 403 (2022)
%         113822 (unconstrained case, P = identity).
%   d_0 = -F_0;  d_k = -F_k + beta_k d_{k-1},  beta_k = c F_k'd_{k-1}/(d_{k-1}'w_{k-1}),
%   w = y + t d_{k-1}, y = F_k - F_{k-1}, t = 1 + max{0, -d_{k-1}'y/||d_{k-1}||^2}.
%   z_k = x_k + alpha_k d_k, alpha_k = max{rho^i}: -F(z_k)'d_k >= sigma*alpha_k*||F(z_k)||*||d_k||^2.
%   x_{k+1} = x_k - xi*gamma_k F(z_k),  gamma_k = F(z_k)'(x_k - z_k)/||F(z_k)||^2.
%   Published parameters: sigma=1e-4, rho=0.55, c=0.6, xi=1.55.
if nargin < 3, opts = struct(); end
[tol, maxiter, maxfev, rec, stopfun] = solver_opts(opts);
def = struct('sigma',1e-4,'rho',0.55,'c',0.6,'xi',1.55);
par = merge_par(def, opts);
tS = tic;
x = x0(:); Fx = F(x); nfev = 1; nbt = 0; k = 0; nF = norm(Fx);
if rec, hist = zeros(maxiter+1,1); hist(1) = nF; end
dp = []; Fprev = []; failed = false;
while nF > tol && k < maxiter && nfev < maxfev
    if k == 0
        d = -Fx;
    else
        yv = Fx - Fprev;  dd = dp'*dp;
        t = 1 + max(0, -(dp'*yv)/dd);
        w = yv + t*dp;
        beta = par.c*(Fx'*dp)/(dp'*w);
        d = -Fx + beta*dp;
    end
    alpha = 1; ok = false; nd2 = d'*d;
    while nfev < maxfev
        z = x + alpha*d;  Fz = F(z);  nfev = nfev + 1;
        if -(Fz'*d) >= par.sigma*alpha*norm(Fz)*nd2
            ok = true; break;
        end
        alpha = par.rho*alpha;  nbt = nbt + 1;
    end
    if ~ok, failed = true; break; end
    if norm(Fz) <= tol
        x = z; Fx = Fz; nF = norm(Fz); k = k + 1;
        if rec, hist(k+1) = nF; end
        break;
    end
    gam = (Fz'*(x - z))/(Fz'*Fz);
    xn = x - par.xi*gam*Fz;
    Fn = F(xn);  nfev = nfev + 1;
    dp = d;  Fprev = Fx;  x = xn;  Fx = Fn;  nF = norm(Fx);  k = k + 1;
    if rec, hist(k+1) = nF; end
    if ~isempty(stopfun) && stopfun(x,k), break; end
end
info.iter = k; info.funcEvals = nfev; info.backtracks = nbt;
info.finalNorm = nF; info.elapsed = toc(tS);
info.success = (nF <= tol); info.failed = failed;
if rec, info.hist = hist(1:k+1); end
end
