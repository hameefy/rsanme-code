function F = logistic_residual(X, c, u, eta)
N = size(X,1);
zu = X*u; r = 1 ./ (1 + exp( c .* zu ));
F = (X' * (-c .* r)) / N + eta * u;
F = F(:);
