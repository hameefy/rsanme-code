function x0 = start_point(id, n)
%START_POINT  The six starting points of the manuscript (test suite).
i = (1:n)';
switch id
    case 1, x0 = 0.8*ones(n,1);
    case 2, x0 = (n - i)/n;            % ((n-1)/n, (n-2)/n, ..., 0)
    case 3, x0 = 0.1*ones(n,1);
    case 4, x0 = i/n;                  % (1/n, 2/n, ..., 1)
    case 5, x0 = 1 - i/n;              % ones - x0^(4)
    case 6, x0 = -ones(n,1);
    otherwise, error('unknown start id');
end
end
