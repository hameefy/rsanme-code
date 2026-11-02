function T = run_monotone(quick)
%RUN_MONOTONE  Main benchmark on Problems 1-5 (dimensions up to 150000, six starting points).
%   Solvers: RSANME and RSANME-Q (tuned parameters from tuned_params.mat if present, otherwise
%   the untuned defaults), Modified Algorithm 2.1, Liu 2022 Algorithm 2.1, DFSANME, PAND-SR
%   (the last four with their published parameters).
%   Budgets (all solvers): ||F(x_k)|| <= 1e-6, at most 1000 iterations and 2000 F evaluations.
%   Output: results_monotone.csv and the performance profiles profile_iter / profile_fev / profile_time.
%   run_monotone(true) = quick smoke test (small dimensions).
if nargin < 1, quick = false; end
root = setup_paths();
outDir = fullfile(root,'results'); if ~exist(outDir,'dir'), mkdir(outDir); end
if quick, dims = [200 500]; else, dims = [1000 5000 10000 50000 100000 150000]; end

opt = struct('tol',1e-6,'maxiter',1000,'maxfev',2000);
parR = struct(); parQ = struct();
f = fullfile(root,'tuned_params.mat');
if exist(f,'file')
    L = load(f);  parR = L.tuned.rsanme;  parQ = L.tuned.rsanmeq;
    fprintf('Using tuned parameters from %s\n', f);
else
    warning('tuned_params.mat not found: RSANME/RSANME-Q run with untuned defaults. Run run_tuning first.');
end
solvers = {@rsanme, @rsanmeq, @modified_algorithm21, @liu2022, @dfsanme, @pandsr};
names   = {'RSANME','RSANME-Q','Modified Alg. 2.1','Liu Alg. 2.1','DFSANME','PAND-SR'};
pars    = {parR, parQ, struct(), struct(), struct(), struct()};

rows = {};  run = 0;
for sid = 1:numel(solvers)
    % one untimed warm-up call so that first-call overhead does not enter the timings
    Fw = test_problem(1);  solvers{sid}(Fw, start_point(1,100), setfield(opt,'par',pars{sid})); %#ok<SFLD>
end
for n = dims
    for pid = 1:5
        [F, pname] = test_problem(pid);
        for gid = 1:6
            x0 = start_point(gid, n);
            for sid = 1:numel(solvers)
                o = opt;  o.par = pars{sid};
                [~, info] = solvers{sid}(F, x0, o);
                run = run + 1;
                rows(end+1,:) = {names{sid}, pname, n, gid, info.iter, info.funcEvals, info.elapsed, info.finalNorm, info.success}; %#ok<AGROW>
            end
            fprintf('%s n=%d start=%d done\n', pname, n, gid);
        end
    end
    T = cell2table(rows, 'VariableNames', {'Solver','Problem','Dim','Start','Iter','Fev','Time','FinalRes','Success'});
    writetable(T, fullfile(outDir,'results_monotone.csv'));      % written after every dimension
end

% ---- performance profiles ------------------------------------------------------------
inst = unique(T(:,{'Problem','Dim','Start'}), 'rows');
metrics = {'Iter','Fev','Time'};  xl = {'iterations','function evaluations','CPU time'};
for m = 1:3
    C = inf(numel(names), height(inst));
    for s = 1:numel(names)
        for i = 1:height(inst)
            sel = strcmp(T.Solver,names{s}) & strcmp(T.Problem,inst.Problem{i}) & T.Dim==inst.Dim(i) & T.Start==inst.Start(i) & T.Success;
            if any(sel), C(s,i) = T.(metrics{m})(find(sel,1)); end
        end
    end
    perfprof(C, names, ['Performance ratio w.r.t. ' xl{m}], fullfile(outDir,['profile_' lower(metrics{m})]));
end
% ---- summary table: solved instances and medians ------------------------------------------
S = table();
for s = 1:numel(names)
    sel = strcmp(T.Solver,names{s});
    S = [S; table(string(names{s}), sum(T.Success(sel)), sum(sel), median(T.Iter(sel & T.Success)), median(T.Fev(sel & T.Success)), median(T.Time(sel & T.Success)), ...
        'VariableNames',{'Solver','Solved','Total','MedianIter','MedianFev','MedianTime'})]; %#ok<AGROW>
end
disp(S);  writetable(S, fullfile(outDir,'summary_monotone.csv'));
end
