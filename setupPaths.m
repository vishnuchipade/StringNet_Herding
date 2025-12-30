% setupPaths.m
% ------------------------------
% Author: Vishnu S. Chipade
% Email: vishnuc@umich.edu
% ------------------------------

function setupPaths()
%SETUPPATHS  Configure MATLAB path for this project

    % Project root = folder containing this file
    projectRoot = fileparts(mfilename('fullpath'));

    % Core code folders (DO add these)
    addpath(fullfile(projectRoot, 'methods'));
    addpath(fullfile(projectRoot, 'config'));
    addpath(genpath(fullfile(projectRoot, 'utils')));
    addpath(fullfile(projectRoot, 'tests'));

    % NOTE:
    % Do NOT add 'data' to the MATLAB path.
    % Access data files via fullfile(projectRoot, 'data', ...)

end


