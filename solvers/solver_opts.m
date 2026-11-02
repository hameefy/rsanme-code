function [tol, maxiter, maxfev, rec, stopfun] = solver_opts(opts)
%SOLVER_OPTS  Common budget / stopping options shared by all solvers.
%   opts.tol      residual tolerance, stop when ||F(x_k)|| <= tol   (default 1e-6)
%   opts.maxiter  iteration budget                                   (default 1000)
%   opts.maxfev   budget on evaluations of F                         (default 2000)
%   opts.record   true -> solver stores the residual history         (default false)
%   opts.stopfun  optional handle  stop = stopfun(x, k)  (extra stopping rule,
%                 used only by the sparse-signal experiment)
if nargin < 1 || isempty(opts), opts = struct(); end
tol = 1e-6;     if isfield(opts,'tol'),     tol = opts.tol;         end
maxiter = 1000; if isfield(opts,'maxiter'), maxiter = opts.maxiter; end
maxfev = 2000;  if isfield(opts,'maxfev'),  maxfev = opts.maxfev;   end
rec = false;    if isfield(opts,'record'),  rec = opts.record;      end
stopfun = [];   if isfield(opts,'stopfun'), stopfun = opts.stopfun; end
end
