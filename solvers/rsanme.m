function [x, info] = rsanme(F, x0, opts)
%RSANME  Positive spectral residual algorithm with a norm-descent line search
%        (Algorithm 1 of the revised manuscript).
%
%   d_0 = -F_0,   d_k = -lambda_k F_k  (k >= 1),
%   lambda_k = tau_k*theta1 + (1-tau_k)*theta2,
%   theta1 = ||s||^2/<y,s>,  theta2 = ||s||/||y||,  y = F_k - F_{k-1} + r*s.
%   Line search: alpha = rho^i, i smallest integer >= 0 with
%       ||F(x_k + alpha d_k)|| <= (1 + eta_k - a*alpha) ||F_k||.
%   alpha is RESET to 1 at every iteration (no warm start).
%
%   Parameters (opts.par overrides):  rho, a, r, eta0, etaq, tau0, tauq
%       eta_k = eta0*etaq^k  (positive, summable  <=>  0<etaq<1, eta0>0)
%       tau_k = tau0*tauq^k  restricted to [1e-12, 1-1e-12]
if nargin < 3, opts = struct(); end
[tol, maxiter, maxfev, rec, stopfun] = solver_opts(opts);
def = struct('rho',0.8,'a',1e-4,'r',1e-3,'eta0',1,'etaq',0.5,'tau0',1,'tauq',0.5);
par = merge_par(def, opts);
tS = tic;
x = x0(:); Fx = F(x); nfev = 1; nbt = 0; k = 0;
nF = norm(Fx);
if rec, hist = zeros(maxiter+1,1); hist(1) = nF; end
xp = []; Fp = [];
failed = false;
while nF > tol && k < maxiter && nfev < maxfev
    if k == 0
        d = -Fx;  lam = 1;
    else
        s = x - xp;  y = Fx - Fp + par.r*s;
        ys = y'*s;  ny = norm(y);  ns = norm(s);
        th1 = ns^2/ys;  th2 = ns/ny;                 % ys >= r*ns^2 > 0 (monotone F)
        tau = min(max(par.tau0*par.tauq^k,1e-12),1-1e-12);
        lam = tau*th1 + (1-tau)*th2;
        d = -lam*Fx;
    end
    eta = par.eta0*par.etaq^k;
    alpha = 1; ok = false;
    while nfev < maxfev
        xt = x + alpha*d;  Ft = F(xt);  nfev = nfev + 1;
        nFt = norm(Ft);
        if nFt <= (1 + eta - par.a*alpha)*nF
            ok = true; break;
        end
        alpha = par.rho*alpha;  nbt = nbt + 1;
    end
    if ~ok, failed = true; break; end
    xp = x;  Fp = Fx;  x = xt;  Fx = Ft;  nF = nFt;  k = k + 1;
    if rec, hist(k+1) = nF; end
    if ~isempty(stopfun) && stopfun(x,k), break; end
end
info.iter = k; info.funcEvals = nfev; info.backtracks = nbt;
info.finalNorm = nF; info.elapsed = toc(tS);
info.success = (nF <= tol); info.failed = failed;
if rec, info.hist = hist(1:k+1); end
end
