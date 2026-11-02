function verify_strong_monotonicity(quick)
%VERIFY_STRONG_MONOTONICITY  For every run of RSANME on Problems 1-5 check, at the computed
%   limit point x*, the hypothesis of the convergence theorem:
%       mu := lambda_min( (J(x*)+J(x*)')/2 ) > 0   (one-sided derivative at kinks),
%   and compare a = (value used) with the threshold mu/(L+r), L = ||J(x*)||_2 estimate
%   (local proxy: ||J(x*)||_2 is used for L).
%   Output: results/strong_monotonicity_check.csv  (mu, ||J||, threshold, a, a < threshold?)
%   The Jacobians are diagonal (P1,P3,P4,P5) or tridiagonal (P2), so lambda_min is computed exactly with eigs/diag.
if nargin < 1, quick = false; end
root = setup_paths();
outDir = fullfile(root,'results'); if ~exist(outDir,'dir'), mkdir(outDir); end
n = 2000; if quick, n = 200; end
par = struct(); f = fullfile(root,'tuned_params.mat');
if exist(f,'file'), L = load(f); par = L.tuned.rsanme; end
if ~isfield(par,'a'), par.a = 1e-4; end
if ~isfield(par,'r'), par.r = 1e-3; end
rows = {};
for pid = 1:5
    [F, pname, Jf] = test_problem(pid);
    for gid = 1:6
        o = struct('par',par,'tol',1e-8);
        [xs, info] = rsanme(F, start_point(gid,n), o);
        J = Jf(xs);  Js = (J + J')/2;
        mu = real(eigs(Js, 1, 'smallestreal', 'Tolerance',1e-10));
        Lest = normest(J);
        thr = mu/(Lest + par.r);
        rows(end+1,:) = {pname, n, gid, info.success, mu, Lest, thr, par.a, par.a < thr}; %#ok<AGROW>
    end
end
T = cell2table(rows,'VariableNames',{'Problem','Dim','Start','Converged','mu','LipschitzEst','mu_over_LplusR','a','a_below_threshold'});
writetable(T, fullfile(outDir,'strong_monotonicity_check.csv'));
disp(T);
fprintf('Hypothesis (mu>0 and a<mu/(L+r)) holds in %d of %d runs.\n', sum(T.mu>0 & T.a_below_threshold), height(T));
end
