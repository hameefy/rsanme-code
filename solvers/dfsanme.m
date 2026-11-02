function [x, info] = dfsanme(F, x0, opts)
%DFSANME  La Cruz, "A spectral algorithm for large-scale systems of nonlinear
%         monotone equations", Numer. Algorithms 76 (2017) 1109-1130, Algorithm 1.
%   x_{k+1} = x_k - lambda_k*alpha_k*F_k, lambda_k = sigma^i (first i >= 0) with
%      ||F(x_k - lambda alpha_k F_k)||^2 <= ||F_k||^2 + eta_k - gamma*lambda^2*||F_k||^2
%   alpha_{k+1} = ||s||^2/<s,y>  if <F_{k+1},F_k>/||F_k||^2 <= 1 - alpha_k/alpha_max,
%               = alpha_max      otherwise;   s = x_{k+1}-x_k, y = F_{k+1}-F_k.
%   Published parameters: gamma=1e-4, sigma=0.5, alpha_0=1, alpha_max=4.5e15,
%   eta_k = 0.999^k (||F_0||^2 + 1e5).
if nargin < 3, opts = struct(); end
[tol, maxiter, maxfev, rec, stopfun] = solver_opts(opts);
def = struct('gamma',1e-4,'sigma',0.5,'alpha0',1,'alphamax',4.5e15,'rhoeta',0.999,'eta0add',1e5);
par = merge_par(def, opts);
tS = tic;
x = x0(:); Fx = F(x); nfev = 1; nbt = 0; k = 0; nF = norm(Fx);
if rec, hist = zeros(maxiter+1,1); hist(1) = nF; end
alk = par.alpha0; failed = false;
eta0 = nF^2 + par.eta0add;
while nF > tol && k < maxiter && nfev < maxfev
    eta = par.rhoeta^k*eta0;
    lam = 1; ok = false;
    while nfev < maxfev
        xt = x - lam*alk*Fx;  Ft = F(xt);  nfev = nfev + 1;
        if norm(Ft)^2 <= nF^2 + eta - par.gamma*lam^2*nF^2
            ok = true; break;
        end
        lam = par.sigma*lam;  nbt = nbt + 1;
    end
    if ~ok, failed = true; break; end
    s = xt - x;  yv = Ft - Fx;  sy = s'*yv;
    if (Ft'*Fx)/nF^2 <= 1 - alk/par.alphamax && sy > 0
        alk = (s'*s)/sy;
    else
        alk = par.alphamax;
    end
    x = xt;  Fx = Ft;  nF = norm(Fx);  k = k + 1;
    if rec, hist(k+1) = nF; end
    if ~isempty(stopfun) && stopfun(x,k), break; end
end
info.iter = k; info.funcEvals = nfev; info.backtracks = nbt;
info.finalNorm = nF; info.elapsed = toc(tS);
info.success = (nF <= tol); info.failed = failed;
if rec, info.hist = hist(1:k+1); end
end
