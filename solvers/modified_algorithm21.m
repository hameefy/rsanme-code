function [x, info] = modified_algorithm21(F, x0, opts)
%MODIFIED_ALGORITHM21  Modified Algorithm 2.1 of K. Muangchoo and A.B. Abubakar,
%   "Enhanced projection method for the solution of the system of nonlinear
%   equations under a more general assumption than pseudo-monotonicity and
%   Lipschitz continuity", Mathematics 12 (2024) 3734, Algorithm 1.
%   Adapted from Experiments/Monotone_Equations/modified_algorithm21.m (same
%   algorithm, same parameters rho=0.8, lambda=1.2, sigma=1e-4); only the
%   interface and the common budget rules were changed.
%   NOTE (inherited from the earlier file): the printed Eq. (12) of the source
%   paper uses F(v_k)'(v_k-u_k); the paper's own proof of Lemma 2 uses
%   F(u_k)'(v_k-u_k).  The proof's version is implemented.
if nargin < 3, opts = struct(); end
[tol, maxiter, maxfev, rec, stopfun] = solver_opts(opts);
def = struct('rho',0.8,'lambda',1.2,'sigma',1e-4);
par = merge_par(def, opts);
tS = tic;
x = x0(:); Fx = F(x); nfev = 1; nbt = 0; k = 0; nF = norm(Fx);
if rec, hist = zeros(maxiter+1,1); hist(1) = nF; end
dp = []; Fprev = []; failed = false;
while nF > tol && k < maxiter && nfev < maxfev
    if k == 0
        d = -Fx;
    else
        yv = Fx - Fprev;  dTy = dp'*yv;  dd = dp'*dp;
        delta = 1 + max(0, -dTy/dd);
        w = yv + delta*dp;
        Fd = Fx'*dp;
        theta = 1 - Fd^2/(nF^2*dd);
        beta = theta*Fd/(dp'*w);
        d = -Fx + beta*dp;
    end
    alpha = 1; ok = false; nd2 = d'*d;
    while nfev < maxfev
        u = x + alpha*d;  Fu = F(u);  nfev = nfev + 1;
        if -(Fu'*d) >= par.sigma*alpha*norm(Fu)*nd2
            ok = true; break;
        end
        alpha = par.rho*alpha;  nbt = nbt + 1;
    end
    if ~ok, failed = true; break; end
    if norm(Fu) <= tol
        x = u; Fx = Fu; nF = norm(Fu); k = k + 1;
        if rec, hist(k+1) = nF; end
        break;
    end
    xi = (Fu'*(x - u))/(Fu'*Fu);
    xn = x - par.lambda*xi*Fu;
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
