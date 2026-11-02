function [x, info] = rsanmeq(F, x0, opts)
%RSANMEQ  Variant of RSANME whose convergence does not restrict the constant a
%         (Algorithm 2 of the revised manuscript).  Same direction as RSANME,
%         quadratic-decrease norm test
%   ||F(x_k + alpha d_k)||^2 <= (1 + eta_k) ||F_k||^2 - sigma*alpha^2*||d_k||^2 .
%   Parameters (opts.par overrides): rho, sigma, r, eta0, etaq, tau0, tauq
if nargin < 3, opts = struct(); end
[tol, maxiter, maxfev, rec, stopfun] = solver_opts(opts);
def = struct('rho',0.8,'sigma',1e-4,'r',1e-3,'eta0',1,'etaq',0.5,'tau0',1,'tauq',0.5);
par = merge_par(def, opts);
tS = tic;
x = x0(:); Fx = F(x); nfev = 1; nbt = 0; k = 0;
nF = norm(Fx);
if rec, hist = zeros(maxiter+1,1); hist(1) = nF; end
xp = []; Fp = [];
failed = false;
while nF > tol && k < maxiter && nfev < maxfev
    if k == 0
        d = -Fx;
    else
        s = x - xp;  y = Fx - Fp + par.r*s;
        ys = y'*s;  ny = norm(y);  ns = norm(s);
        th1 = ns^2/ys;  th2 = ns/ny;
        tau = min(max(par.tau0*par.tauq^k,1e-12),1-1e-12);
        d = -(tau*th1 + (1-tau)*th2)*Fx;
    end
    eta = par.eta0*par.etaq^k;  nd2 = d'*d;
    alpha = 1; ok = false;
    while nfev < maxfev
        xt = x + alpha*d;  Ft = F(xt);  nfev = nfev + 1;
        nFt = norm(Ft);
        if nFt^2 <= (1 + eta)*nF^2 - par.sigma*alpha^2*nd2
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
