function root = setup_paths()
%SETUP_PATHS  Add solvers/, problems/ and drivers/ to the MATLAB path; returns the code root.
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root,'solvers'), fullfile(root,'problems'), fullfile(root,'drivers'));
end
