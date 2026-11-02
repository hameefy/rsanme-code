function par = merge_par(defaults, opts)
%MERGE_PAR  Overwrite the fields of DEFAULTS by those of opts.par (if present).
par = defaults;
if nargin >= 2 && isstruct(opts) && isfield(opts,'par') && isstruct(opts.par)
    f = fieldnames(opts.par);
    for i = 1:numel(f), par.(f{i}) = opts.par.(f{i}); end
end
end
