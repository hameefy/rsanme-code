function [x, info] = pandsr(F, x0, opts)
%PANDSR  Pand-SR of Morini, Porcelli and Toint, "Approximate norm descent
%        methods for constrained nonlinear systems", Math. Comp. 87 (2018)
%        1327-1351 (Algorithm 3.1 with B_k = beta_k^{-1} I, i.e. the spectral
%        residual step of Psane), specialised to Omega = R^n (P = identity).
%   p^QN = -beta_k F_k;  trial steps p_+ = lam*p^QN,  p_- = -lam*p^QN.
%   (14)  ||F(x+p)|| <= (1 - alpha(1+lam)) ||F_k||       tested for p_+, then p_-
%   (15)  ||F(x+p)|| <= (1 + eta_k - alpha*lam) ||F_k||  tested for p_+, then p_-
%         (nonzero steps only);  otherwise lam <- sigma*lam.
%   beta_{k+1} = 1/b_k if |1/b_k| in [bmin,bmax], else |1/b_k| projected onto
%   [bmin,bmax];   b_k = s'y/s's.
%   Published parameters: bmin=1e-30, bmax=1e30, beta_0=1, alpha=1e-4,
%   sigma=0.5, lam_max=1, eta_k = 0.99^k (100 + ||F_0||^2).
if nargin < 3, opts = struct(); end
[tol, maxiter, maxfev, rec, stopfun] = solver_opts(opts);
def = struct('alpha',1e-4,'sigma',0.5,'bmin',1e-30,'bmax',1e30,'beta0',1,'rhoeta',0.99,'eta0add',100);
par = merge_par(def, opts);
tS = tic;
x = x0(:); Fx = F(x); nfev = 1; nbt = 0; k = 0; nF = norm(Fx);
if rec, hist = zeros(maxiter+1,1); hist(1) = nF; end
beta = par.beta0; failed = false;
eta0 = nF^2 + par.eta0add;
while nF > tol && k < maxiter && nfev < maxfev
    eta = par.rhoeta^k*eta0;
    pqn = -beta*Fx;
    lam = 1; ok = false;
    while nfev < maxfev
        pp = lam*pqn;  pm = -lam*pqn;
        xp_ = x + pp;  Fp_ = F(xp_);  nfev = nfev + 1;
        if norm(Fp_) <= (1 - par.alpha*(1+lam))*nF
            xt = xp_; Ft = Fp_; ok = true; break;
        end
        if nfev >= maxfev, break; end
        xm_ = x + pm;  Fm_ = F(xm_);  nfev = nfev + 1;
        if norm(Fm_) <= (1 - par.alpha*(1+lam))*nF
            xt = xm_; Ft = Fm_; ok = true; break;
        end
        if norm(pp) > 0 && norm(Fp_) <= (1 + eta - par.alpha*lam)*nF
            xt = xp_; Ft = Fp_; ok = true; break;
        end
        if norm(pm) > 0 && norm(Fm_) <= (1 + eta - par.alpha*lam)*nF
            xt = xm_; Ft = Fm_; ok = true; break;
        end
        lam = par.sigma*lam;  nbt = nbt + 1;
    end
    if ~ok, failed = true; break; end
    s = xt - x;  yv = Ft - Fx;  bk = (s'*yv)/(s'*s);
    ib = 1/bk;
    if abs(ib) >= par.bmin && abs(ib) <= par.bmax
        beta = ib;
    else
        beta = min(par.bmax, max(par.bmin, abs(ib)));
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
