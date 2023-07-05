% Get the full path to the root directory
root_directory = pwd; % Assuming you are executing startup.m from Inverse-Nonlinear-Schrodinger

% Specify the directories within the root_directory
src_directory = fullfile(root_directory, 'src');
test_directory = fullfile(root_directory, 'tests');
examples_directory = fullfile(root_directory, 'examples');

% Add the directories to the MATLAB search path
addpath(src_directory);
addpath(test_directory);
addpath(examples_directory);
