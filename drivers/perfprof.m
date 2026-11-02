function rho = perfprof(C, names, xlab, outfile, tauMax)
%PERFPROF  Dolan-More performance profile.
%   C      : nSolvers x nProblems cost matrix; failures must be Inf (or NaN)
%   names  : cellstr of solver names
%   xlab   : label of the horizontal axis
%   outfile: base name; writes <outfile>.png and <outfile>.eps (skipped if empty)
%   rho(s,t) = fraction of problems with  C(s,p) <= tau_t * min_s C(s,p)
if nargin < 5, tauMax = 1e3; end
C(~isfinite(C)) = Inf;
best = min(C,[],1);
keep = isfinite(best);                       % problems solved by at least one solver
C = C(:,keep); best = best(keep);
taus = unique([1, logspace(0, log10(tauMax), 400)]);
rho = zeros(size(C,1), numel(taus));
for t = 1:numel(taus)
    rho(:,t) = mean(C <= taus(t)*best, 2);
end
h = figure('Color','w','Visible','off');
sty = {'-','--','-.',':','-','--'};
hold on;
for s = 1:size(C,1)
    stairs(taus, rho(s,:), sty{mod(s-1,numel(sty))+1}, 'LineWidth', 1.8);
end
set(gca,'XScale','log'); grid on; box on;
xlabel(xlab); ylabel('P(\tau)');
xlim([1 tauMax]); ylim([0 1.02]);
legend(names, 'Location','southeast');
if ~isempty(outfile)
    exportgraphics(h, [outfile '.png'], 'Resolution', 300);
    print(h, [outfile '.eps'], '-depsc');
end
close(h);
end
