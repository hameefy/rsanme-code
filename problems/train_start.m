function x0 = train_start(set, id, n)
%TRAIN_START  Starting points for tuning ('train') and validation ('valid') -- disjoint from the test set.
i = (1:n)';
if strcmp(set,'train')
    switch id
        case 1, x0 = 2*ones(n,1);
        case 2, x0 = -0.5*ones(n,1);
        case 3, x0 = linspace(-1,2,n)';
        case 4, s = RandStream('mt19937ar','Seed',7);  x0 = 3*(rand(s,n,1) - 0.3);
    end
else
    switch id
        case 1, x0 = 3*ones(n,1);
        case 2, x0 = -2*ones(n,1);
        case 3, x0 = 2*cos(i);
        case 4, s = RandStream('mt19937ar','Seed',11);  x0 = 4*(rand(s,n,1) - 0.5);
    end
end
end
