function run_tuning(quick)
%RUN_TUNING  Tune the free parameters of RSANME and RSANME-Q on a TRAINING suite that is
%   disjoint from the test suite (different problems TR1-TR5, dimensions and starting points),
%   check the winner on a separate VALIDATION suite, and write tuned_params.mat / .txt.
%
%   Score of a configuration (lower is better), computed over all training instances i:
%       score = mean_i log2( min( cost_i / best_i , 100 ) ),     failure -> ratio 100,
%   cost = number of F evaluations, best_i = best cost among ALL configurations tried on i.
%   Failures (||F|| > 1e-6 within 1000 iterations / 2000 F evaluations) are therefore
%   penalised first and efficiency second.  Stage 1: random search (log-uniform);
%   Stage 2: coordinate refinement around the 5 best.  Final choice: smallest
%   VALIDATION score among the 5 refined candidates (guards against over-fitting).
%
%   run_tuning(true) runs a very small smoke version.
if nargin < 1, quick = false; end
root = setup_paths();  rng(2026,'twister');
if quick, dimsT = 200; dimsV = 300; nRand = 6; nRefine = 1; else, dimsT = [1000 10000]; dimsV = [2000 20000]; nRand = 300; nRefine = 2; end

algs = struct('name',{'rsanme','rsanmeq'}, 'fh',{@rsanme,@rsanmeq});
for ai = 1:numel(algs)
    name = algs(ai).name;  fh = algs(ai).fh;
    fprintf('=== Tuning %s ===\n', name);
    trainSet = build_set('train', dimsT);   validSet = build_set('valid', dimsV);
    space = param_space(name);
    % ---- Stage 1: random search -------------------------------------------------
    cfgs = cell(nRand+1,1);  cfgs{1} = default_cfg(name);
    for c = 2:nRand+1, cfgs{c} = sample_cfg(space); end
    COST = eval_cfgs(fh, cfgs, trainSet);           % nCfg x nInst
    bestVec = min(COST,[],1);
    sc = scores(COST, bestVec);
    [~,ord] = sort(sc);  top = cfgs(ord(1:5))';
    fprintf('stage 1 best score %.3f (default %.3f)\n', sc(ord(1)), sc(1));
    % ---- Stage 2: coordinate refinement (training set only) ---------------------
    fn = fieldnames(space);
    for t = 1:5
        cur = top{t};
        for rep = 1:nRefine
            for p = 1:numel(fn)
                f = space.(fn{p});  cand = cell(1,5);
                for j = 1:5
                    cand{j} = cur;  val = cur.(fn{p});
                    if f.log, newv = val*10^((j-3)*0.25); else, newv = val + (j-3)*0.05*(f.hi-f.lo); end
                    cand{j}.(fn{p}) = min(max(newv, f.lo), f.hi);
                end
                M = eval_cfgs(fh, cand, trainSet);
                bestVec = min([bestVec; M],[],1);
                [~,bj] = min(scores(M, bestVec));
                cur = cand{bj};
            end
        end
        top{t} = cur;
    end
    % ---- Validation --------------------------------------------------------------
    Vcost = eval_cfgs(fh, [top, {default_cfg(name)}], validSet);
    vs = scores(Vcost, min(Vcost,[],1));
    [~,bv] = min(vs(1:5));
    best = top{bv};
    Tcost = eval_cfgs(fh, [top, {default_cfg(name)}], trainSet);  ts = scores(Tcost, min(Tcost,[],1));
    fprintf('validation score: tuned %.3f  default %.3f\n', vs(bv), vs(6));
    fprintf('training  score: tuned %.3f  default %.3f\n', ts(bv), ts(6));
    fprintf('failures (valid): tuned %d  default %d  of %d\n', sum(~isfinite(Vcost(bv,:))), sum(~isfinite(Vcost(6,:))), size(Vcost,2));
    tuned.(name) = best;  report.(name) = struct('train_tuned',ts(bv),'train_default',ts(6),'valid_tuned',vs(bv),'valid_default',vs(6));
    disp(best);
end
save(fullfile(root,'tuned_params.mat'), 'tuned', 'report');
fid = fopen(fullfile(root,'tuned_params.txt'),'w');
for ai = 1:numel(algs)
    n = algs(ai).name; fprintf(fid,'%s\n', n);
    f = fieldnames(tuned.(n));
    for i = 1:numel(f), fprintf(fid,'  %-8s = %.6g\n', f{i}, tuned.(n).(f{i})); end
    r = report.(n);
    fprintf(fid,'  score train: tuned %.3f default %.3f | valid: tuned %.3f default %.3f\n', r.train_tuned, r.train_default, r.valid_tuned, r.valid_default);
end
fclose(fid);
fprintf('Wrote tuned_params.mat and tuned_params.txt in %s\n', root);
end

% ---------------------------------------------------------------------------------
function space = param_space(name)
% All admissible by the theory: rho in (0,1), a in (0,1), r>0, eta0>0, 0<etaq<1, tau in (0,1).
space.rho  = struct('lo',0.3,  'hi',0.95, 'log',false);
space.r    = struct('lo',1e-5, 'hi',1,    'log',true);
space.eta0 = struct('lo',1e-3, 'hi',1e3,  'log',true);
space.etaq = struct('lo',0.3,  'hi',0.99, 'log',false);
space.tau0 = struct('lo',0.05, 'hi',0.95, 'log',false);
space.tauq = struct('lo',0.3,  'hi',1.0,  'log',false);
if strcmp(name,'rsanme')
    space.a = struct('lo',1e-7,'hi',0.5,'log',true);
else
    space.sigma = struct('lo',1e-7,'hi',1,'log',true);
end
end

function c = default_cfg(name)
c = struct('rho',0.8,'r',1e-3,'eta0',1,'etaq',0.5,'tau0',1,'tauq',0.5);
c.tau0 = 0.95;                                   % tau_k = 1/2^k of the manuscript is tau0=1 -> clipped to <1
if strcmp(name,'rsanme'), c.a = 1e-4; else, c.sigma = 1e-4; end
end

function c = sample_cfg(space)
fn = fieldnames(space);  c = struct();
for i = 1:numel(fn)
    f = space.(fn{i});
    if f.log, c.(fn{i}) = 10^(log10(f.lo) + rand*(log10(f.hi)-log10(f.lo)));
    else,     c.(fn{i}) = f.lo + rand*(f.hi - f.lo); end
end
end

function set = build_set(which, dims)
set = struct('F',{},'x0',{});  k = 0;
for n = dims
    for pid = 1:5
        if pid == 3 && n > 2000, continue; end          % dense affine problem
        F = train_problem(pid, n);
        for sid = 1:4
            k = k + 1;  set(k).F = F;  set(k).x0 = train_start(which, sid, n);
        end
    end
end
end

function COST = eval_cfgs(fh, cfgs, set)
COST = inf(numel(cfgs), numel(set));
for c = 1:numel(cfgs)
    for i = 1:numel(set)
        [~,info] = fh(set(i).F, set(i).x0, struct('par',cfgs{c}));
        if info.success, COST(c,i) = info.funcEvals; end
    end
end
end

function s = scores(COST, best)
R = COST ./ best;  R(~isfinite(R)) = 100;  R(isnan(R)) = 100;  R = min(R,100);
s = mean(log2(R),2);
end
