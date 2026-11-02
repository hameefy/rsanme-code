function run_logistic(quick)
%RUN_LOGISTIC  l2-regularised logistic regression on LIBSVM data (Section "Logistic regression").
%   F(u) = (1/Z) sum_i -m_i exp(-m_i n_i'u) n_i / (1 + exp(-m_i n_i'u)) + gamma*u,  gamma = 0.1.
%   F is gamma-strongly monotone (gradient of a gamma-strongly convex function) and Lipschitz, so the
%   hypotheses of the convergence theorem hold with mu = gamma, globally.
%   Solvers: RSANME, RSANME-Q, Modified Algorithm 2.1, ISCGPM, DFSANME, PAND-SR.
%   12 random starting points u0 = 5*(rand(n,1)-0.1) per data set (same u0 for all solvers).
%   Datasets are read from DATADIR (edit below); libsvmread.mex* must be on the path (it is shipped in
%   Experiments/Logistic_Regression).
if nargin < 1, quick = false; end
root = setup_paths();
dataDir = fullfile(root,'..','..','Experiments','Logistic_Regression');   % <-- edit if the data live elsewhere
addpath(dataDir);
outDir = fullfile(root,'results'); if ~exist(outDir,'dir'), mkdir(outDir); end
files = {'w1a.t','w2a.t','w3a.t','w4a.t','w5a.t','w6a.t','w7a.t','w8a.t','svmguide4.t','phishing.txt','fourclass_scale.txt','a9a.t'};
if quick, files = files([1 9 11]); end
reps = 12; if quick, reps = 2; end
gam = 0.1;  rng(42,'twister');
par = struct(); parq = struct();
f = fullfile(root,'tuned_params.mat');
if exist(f,'file'), L = load(f); par = L.tuned.rsanme; parq = L.tuned.rsanmeq; end
names = {'RSANME','RSANME-Q','Modified Alg. 2.1','ISCGPM','DFSANME','PAND-SR'};
opt = struct('tol',1e-6,'maxiter',1000,'maxfev',2000,'record',true);
rows = {};  hist = cell(numel(files), numel(names));
for di = 1:numel(files)
    [Y, X] = libsvmread(fullfile(dataDir, files{di}));
    u = unique(Y);
    if ~all(ismember(u,[-1 1])), s = sort(u); Y(Y==s(1)) = -1; Y(Y==s(end)) = 1; end
    X = sparse(X);  n = size(X,2);  Z = size(X,1);
    Fh = @(w) logistic_residual(X, Y, w, gam);
    fprintf('--- %s (Z=%d, n=%d)\n', files{di}, Z, n);
    for r = 1:reps
        u0 = 5*(rand(n,1) - 0.1);
        for s = 1:numel(names)
            switch names{s}
                case 'RSANME',           o = opt; o.par = par;  [~,info] = rsanme(Fh,u0,o);
                case 'RSANME-Q',         o = opt; o.par = parq; [~,info] = rsanmeq(Fh,u0,o);
                case 'Modified Alg. 2.1', [~,info] = modified_algorithm21(Fh,u0,opt);
                case 'DFSANME',          [~,info] = dfsanme(Fh,u0,opt);
                case 'PAND-SR',          [~,info] = pandsr(Fh,u0,opt);
                case 'ISCGPM'
                    p = struct('alphaBar',0.7,'etaSpec',0.3,'nu',0.785,'mu',6,'r',0.26,'gamma',1.9,'rho',0.32,'sigma',1e-4,'tol',1e-6,'maxIter',1000);
                    t0 = tic; [xI,ii] = ISCGPM_solver(u0, X, Y, gam, p); el = toc(t0);
                    fin = norm(Fh(xI));                  % residual at the RETURNED point
                    info = struct('iter',ii.iter,'funcEvals',ii.fval,'elapsed',el,'finalNorm',fin,'success',(ii.iter<1000),'hist',ii.resHist);
            end
            rows(end+1,:) = {files{di}, names{s}, r, info.iter, info.funcEvals, info.elapsed, info.finalNorm, info.success}; %#ok<AGROW>
            if r == 1, hist{di,s} = info.hist; end
        end
    end
    T = cell2table(rows,'VariableNames',{'Dataset','Solver','Rep','Iter','Fev','Time','FinalRes','Success'});
    writetable(T, fullfile(outDir,'results_logistic_runs.csv'));
end
% mean over repeats (as in the manuscript table)
G = groupsummary(T, {'Dataset','Solver'}, 'mean', {'Iter','Fev','Time','FinalRes','Success'});
writetable(G, fullfile(outDir,'results_logistic_summary.csv'));  disp(G);
% residual histories (first repetition) : relative residual vs iteration
for di = 1:numel(files)
    h = figure('Color','w','Visible','off'); hold on;
    for s = 1:numel(names)
        v = hist{di,s}; if isempty(v), continue; end
        semilogy(0:numel(v)-1, v/v(1), 'LineWidth', 1.8);
    end
    set(gca,'YScale','log'); grid on; box on; xlabel('iteration'); ylabel('||F(u_k)||/||F(u_0)||');
    legend(names,'Location','northeast'); title(strrep(files{di},'_','\_'));
    [~,base] = fileparts(files{di});
    exportgraphics(h, fullfile(outDir,['Residual_' base '.png']), 'Resolution', 300);
    print(h, fullfile(outDir,['Residual_' base '.eps']), '-depsc');  close(h);
end
end
