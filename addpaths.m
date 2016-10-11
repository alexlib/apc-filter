function addpaths(ROOT)
% Root is the top level directory containing all
% the other code directories


% Get the path to the user's desktop
if nargin < 1
    ROOT = fullfile(getenv('HOME'), 'Desktop');
end

% Add the paths
addpath(fullfile(ROOT, 'spectral-phase-correlation', 'filtering'));
addpath(fullfile(ROOT, 'spectral-phase-correlation', 'correlation_algorithms'))
addpath(fullfile(ROOT, 'piv-image-generation'));
addpath(fullfile(ROOT, 'piv-image-generation', 'jobfiles'));

% Add the other libraries if you're on mac
if ismac
    addpath(fullfile(ROOT, 'tightfig'));
    addpath(fullfile(ROOT, 'shadedErrorBar'));
    addpath(fullfile(ROOT, 'export_fig'));
end

end