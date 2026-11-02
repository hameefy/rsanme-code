function [Fh, name] = train_problem(id, n)
%TRAIN_PROBLEM  Monotone Lipschitz problems used ONLY for tuning/validating RSANME
%   parameters.  They are disjoint from the test problems P1-P5.
%   TR1  F_i = 2 x_i + sin x_i                         (smooth, strongly monotone)
%   TR2  F_i = x_i - sin|x_i|                          (nonsmooth, degenerate root: NOT strongly monotone)
%   TR3  F(x) = (N N' + U + D) x / n   (random affine, U skew)  [n <= 2000]
%   TR4  tridiagonal  F_i = 3x_i - x_{i-1} - x_{i+1} + 0.5(x_{i+1}-x_{i-1})  (skew part added)
%   TR5  F_i = 0.1 x_i + max(x_i,0)                    (nonsmooth, strongly monotone)
switch id
    case 1
        name = 'TR1';  Fh = @(x) 2*x + sin(x);
    case 2
        name = 'TR2';  Fh = @(x) x - sin(abs(x));
    case 3
        name = 'TR3';
        s = RandStream('mt19937ar','Seed',100+n);
        N = 4*rand(s,n,n) - 2;  U1 = 4*rand(s,n,n) - 2;
        U = 0.5*(U1' - U1);  D = diag(2*rand(s,n,1));
        M = (N*N' + U + D)/n;   %#ok<MINV>
        Fh = @(x) M*x;
    case 4
        name = 'TR4';
        Fh = @tri;
    case 5
        name = 'TR5';  Fh = @(x) 0.1*x + max(x,0);
    otherwise, error('unknown training problem');
end
end

function F = tri(x)
xm = [0; x(1:end-1)];  xp = [x(2:end); 0];
F = 3*x - xm - xp + 0.5*(xp - xm);
end
