function [Fh, name, Jdiag] = test_problem(id)
%TEST_PROBLEM  Benchmark Problems 1-5 of the manuscript (Section "Monotone Nonlinear Equations").
%   Fh    : handle F(x) -> column vector
%   Jdiag : handle returning the (generalised) Jacobian as a sparse matrix, used only by
%           verify_strong_monotonicity.m (one-sided derivative where F is not differentiable).
switch id
    case 1   % Yu et al. 2009:  F_i = x_i - sin|x_i - 1|
        name = 'P1';
        Fh = @(x) x - sin(abs(x-1));
        Jdiag = @(x) spdiags(1 - cos(abs(x-1)).*sign(x-1), 0, numel(x), numel(x));
    case 2   % Yu et al. 2013: tridiagonal exponential-cosine
        name = 'P2';
        Fh = @p2;
        Jdiag = @p2jac;
    case 3   % Zhou-Li 2007:  F_i = 2x_i - sin|x_i|
        name = 'P3';
        Fh = @(x) 2*x - sin(abs(x));
        Jdiag = @(x) spdiags(2 - cos(x).*sign(x), 0, numel(x), numel(x));
    case 4   % La Cruz 2014:  F_i = cos(x_i) + x_i - 1
        name = 'P4';
        Fh = @(x) cos(x) + x - 1;
        Jdiag = @(x) spdiags(1 - sin(x), 0, numel(x), numel(x));
    case 5   % Wang et al. 2007:  F_i = exp(x_i) - 1
        name = 'P5';
        Fh = @(x) exp(x) - 1;
        Jdiag = @(x) spdiags(exp(x), 0, numel(x), numel(x));
    otherwise
        error('unknown problem id');
end
end

function F = p2(x)
n = numel(x);  h = 1/(n+1);
t = zeros(n,1);
t(1)     = x(1) + x(2);
t(n)     = x(n-1) + x(n);
t(2:n-1) = x(1:n-2) + x(2:n-1) + x(3:n);
F = x - exp(cos(h*t));
end

function J = p2jac(x)
n = numel(x);  h = 1/(n+1);
t = zeros(n,1);
t(1)     = x(1) + x(2);
t(n)     = x(n-1) + x(n);
t(2:n-1) = x(1:n-2) + x(2:n-1) + x(3:n);
g = h*sin(h*t).*exp(cos(h*t));     % dF_i/dx_j = delta_ij + g_i for j in the stencil of row i
J = speye(n);
for d = -1:1
    idx = (max(1,1-d):min(n,n-d))';
    J = J + sparse(idx, idx+d, g(idx), n, n);
end
end
