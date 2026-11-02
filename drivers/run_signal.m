function run_signal(quick)
%RUN_SIGNAL  Sparse signal reconstruction, min 0.5||Du-b||^2 + gamma||u||_1, via the monotone map
%   F(u) = min{u, Bu + c}  (u = [u+; u-]), B = [D'D -D'D; -D'D D'D], c = gamma*e +/- D'b.
%   F is monotone and Lipschitz but NOT strongly monotone (B is singular), so this experiment is
%   outside the hypotheses of the convergence theorem; it is reported as a practical illustration.
%   Setting taken from the manuscript text: D = randn(k,n), noise variance 1e-12, gamma = 0.5,
%   u0 = D'b (split), stopping rule ||F(u_k)|| <= 1e-6 ||F(u_0)|| (the relative-objective rule of the first version stopped all
%   runs after ~25 iterations with large residuals and was dropped), budgets 5000 iterations / 10000 F-evaluations.
%   n = 4096, k = 1024, 64 non-zeros, 10 random instances (confirmed by the author).
%   Solvers: RSANME, RSANME-Q, Modified Alg. 2.1, DFSANME, PAND-SR.  Metrics: iterations, F evaluations,
%   time, MSE, SNR, residual ||F||.
if nargin < 1, quick = false; end
root = setup_paths();
outDir = fullfile(root,'results'); if ~exist(outDir,'dir'), mkdir(outDir); end
n = 4096; k = 1024; nnzs = 64; ninst = 10; gam = 0.5; sig2 = 1e-12;
if quick, n = 512; k = 128; nnzs = 8; ninst = 2; end
par = struct(); parq = struct();
f = fullfile(root,'tuned_params.mat');
if exist(f,'file'), L = load(f); par = L.tuned.rsanme; parq = L.tuned.rsanmeq; end
names = {'RSANME','RSANME-Q','Modified Alg. 2.1','DFSANME','PAND-SR'};
rows = {};  rng(2026,'twister');
for inst = 1:ninst
    ub = zeros(n,1);  idx = randperm(n, nnzs);  ub(idx) = randn(nnzs,1);
    D = randn(k,n);  b = D*ub + sqrt(sig2)*randn(k,1);
    DtD = D'*D;  Dtb = D'*b;
    cvec = [gam + (-Dtb); gam + Dtb];
    Bmul = @(z) [DtD*(z(1:n) - z(n+1:end)); -DtD*(z(1:n) - z(n+1:end))];
    Fh = @(z) min(z, Bmul(z) + cvec);
    obj = @(z) 0.5*norm(D*(z(1:n)-z(n+1:end)) - b)^2 + gam*sum(z);
    u0 = [max(Dtb,0); max(-Dtb,0)];                         % split of D'b
    u0n = norm(Fh(u0));
    o = struct('tol',1e-6*u0n,'maxiter',5000,'maxfev',10000);   % residual rule ||F|| <= 1e-6 ||F(u0)||, larger budget
    for s = 1:numel(names)
        oo = o;
        switch names{s}
            case 'RSANME',           oo.par = par;  [z,info] = rsanme(Fh,u0,oo);
            case 'RSANME-Q',         oo.par = parq; [z,info] = rsanmeq(Fh,u0,oo);
            case 'Modified Alg. 2.1', [z,info] = modified_algorithm21(Fh,u0,oo);
            case 'DFSANME',          [z,info] = dfsanme(Fh,u0,oo);
            case 'PAND-SR',          [z,info] = pandsr(Fh,u0,oo);
        end
        rec = z(1:n) - z(n+1:end);
        mse = norm(ub - rec)^2/n;
        snr = 10*log10(norm(ub)^2/norm(ub - rec)^2);
        rows(end+1,:) = {inst, names{s}, info.iter, info.funcEvals, info.elapsed, norm(Fh(z)), mse, snr, obj(z), double(norm(Fh(z)) <= o.tol)}; %#ok<AGROW>
        if inst == 1, save(fullfile(outDir,sprintf('signal_%s.mat',matlab.lang.makeValidName(names{s}))),'rec','ub'); end
    end
    fprintf('instance %d done\n', inst);
end
T = cell2table(rows,'VariableNames',{'Instance','Solver','Iter','Fev','Time','ResF','MSE','SNR_dB','Objective','Success'});
writetable(T, fullfile(outDir,'results_signal_runs.csv'));
G = groupsummary(T,'Solver','mean',{'Iter','Fev','Time','ResF','MSE','SNR_dB','Objective','Success'});
writetable(G, fullfile(outDir,'results_signal_summary.csv'));  disp(G);
end

function sf = make_stop(obj, f0, tol)
% returns a handle with memory of the previous objective value (persistent via nested closure)
store = struct('prev', f0);
sf = @(z,kk) stopstep(z, kk);
    function st = stopstep(z, ~)
        fz = obj(z);
        st = abs(fz - store.prev) < tol*abs(store.prev);
        store.prev = fz;
    end
end
