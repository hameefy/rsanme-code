function run_sensitivity_a(quick)
%RUN_SENSITIVITY_A  Effect of the line-search constant a on RSANME (the theory requires
%   a < mu/(L+r) at the limit point) and, for comparison, of sigma on RSANME-Q (no restriction).
%   Problems 1-5, n = 10000, six starting points; all other RSANME parameters at their tuned values.
%   Reports per value: solved instances, mean iterations, mean F evaluations, mean accepted step
%   size and the fraction of accepted steps with alpha_k < 1e-3.
if nargin < 1, quick = false; end
root = setup_paths();
outDir = fullfile(root,'results'); if ~exist(outDir,'dir'), mkdir(outDir); end
n = 10000; if quick, n = 300; end
par = struct();  parq = struct();
f = fullfile(root,'tuned_params.mat');
if exist(f,'file'), L = load(f); par = L.tuned.rsanme; parq = L.tuned.rsanmeq; end
avals = [1e-8 1e-6 1e-4 1e-3 1e-2 1e-1 0.3 0.6 0.9];
opt = struct('tol',1e-6,'maxiter',1000,'maxfev',2000);
rows = {};
for a = avals
    solved = 0; its = []; fevs = []; amean = []; small = [];
    for pid = 1:5
        F = test_problem(pid);
        for gid = 1:6
            p = par; p.a = a; o = opt; o.par = p;
            [~, info] = rsanme_trace(F, start_point(gid,n), o);
            solved = solved + info.success;  its(end+1) = info.iter; fevs(end+1) = info.funcEvals; %#ok<AGROW>
            amean(end+1) = info.meanAlpha; small(end+1) = info.fracSmall; %#ok<AGROW>
        end
    end
    rows(end+1,:) = {'RSANME','a',a,solved,30,mean(its),mean(fevs),mean(amean),mean(small)}; %#ok<AGROW>
    fprintf('a=%g solved %d/30 mean it %.1f fev %.1f mean alpha %.3f\n', a, solved, mean(its), mean(fevs), mean(amean));
end
sigvals = [1e-8 1e-6 1e-4 1e-2 1];
for sg = sigvals
    solved = 0; its = []; fevs = [];
    for pid = 1:5
        F = test_problem(pid);
        for gid = 1:6
            p = parq; p.sigma = sg; o = opt; o.par = p;
            [~, info] = rsanmeq(F, start_point(gid,n), o);
            solved = solved + info.success; its(end+1) = info.iter; fevs(end+1) = info.funcEvals; %#ok<AGROW>
        end
    end
    rows(end+1,:) = {'RSANME-Q','sigma',sg,solved,30,mean(its),mean(fevs),NaN,NaN}; %#ok<AGROW>
    fprintf('sigma=%g solved %d/30 mean it %.1f fev %.1f\n', sg, solved, mean(its), mean(fevs));
end
T = cell2table(rows,'VariableNames',{'Solver','Param','Value','Solved','Total','MeanIter','MeanFev','MeanAlpha','FracAlphaBelow1e3'});
writetable(T, fullfile(outDir,'sensitivity_a.csv')); disp(T);
end

function [x, info] = rsanme_trace(F, x0, opts)
% RSANME with the accepted step sizes recorded (copy of rsanme.m + trace of alpha_k).
[tol, maxiter, maxfev] = solver_opts(opts);
def = struct('rho',0.8,'a',1e-4,'r',1e-3,'eta0',1,'etaq',0.5,'tau0',1,'tauq',0.5);
par = merge_par(def, opts);
x = x0(:); Fx = F(x); nfev = 1; k = 0; nF = norm(Fx); xp = []; Fp = []; al = zeros(maxiter,1);
while nF > tol && k < maxiter && nfev < maxfev
    if k == 0, d = -Fx;
    else
        s = x - xp; y = Fx - Fp + par.r*s; ns = norm(s);
        tau = min(max(par.tau0*par.tauq^k,1e-12),1-1e-12);
        d = -(tau*ns^2/(y'*s) + (1-tau)*ns/norm(y))*Fx;
    end
    eta = par.eta0*par.etaq^k; alpha = 1; ok = false;
    while nfev < maxfev
        xt = x + alpha*d; Ft = F(xt); nfev = nfev + 1;
        if norm(Ft) <= (1 + eta - par.a*alpha)*nF, ok = true; break; end
        alpha = par.rho*alpha;
    end
    if ~ok, break; end
    xp = x; Fp = Fx; x = xt; Fx = Ft; nF = norm(Fx); k = k + 1; al(k) = alpha;
end
info.iter = k; info.funcEvals = nfev; info.finalNorm = nF; info.success = nF <= tol;
if k > 0, info.meanAlpha = mean(al(1:k)); info.fracSmall = mean(al(1:k) < 1e-3);
else, info.meanAlpha = NaN; info.fracSmall = NaN; end
end
