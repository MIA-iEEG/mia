% ========================================================================
% This file is part of MIA.
% 
% MIA is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% MIA is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%  
% Copyright (C) 2016-2025 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)
% 


%% TODO: Automate the creation of the data structure (dat)
% % Input files
% sFiles = [];
% SubjectNames = {...
%     'COREG'};
% 
% chan = load('/Users/annesophiedubarry/Library/CloudStorage/SynologyDrive-NAS/0_projects/in_progress/EXPRESSONINTRA/data/EXPRESSONINTRA/data/COREG/@intra/channel.mat');
% 
% % Get the total number of channels
% Nchan = length(chan.Channel); 
% 
% % Process: Simulate generic signals
% sFiles = bst_process('CallProcess', 'process_simulate_matrix', sFiles, [], ...
%     'subjectname', SubjectNames{1}, ...
%     'condition',   'Simulation', ...
%     'samples',     10, ...
%     'srate',       1000, ...
%     'matlab',      ['Data(1,:) = sin(2*pi*t);' 10 strcat('Data(',num2str(Nchan),',:) = cos(pi*t) + 1;')]);

% ========================================================================

% Before executing this code, export channels and data from COREG.
% To do this, right-click on Channel and Data > File -> Export to Matlab.

%% Color Code by Patient
% Extract a vector of patient IDs (digits following the letter "S")
vector = cellfun(@(x) sscanf(x, 'S%d'), {chan.Channel.Comment}, 'UniformOutput', true);

% Get unique patient IDs and map them to sequential integers
[unique_vals, ~, idx] = unique(vector);

% Replace all data in dat.F with a patient-based color code
dat.F = repmat(idx, 1, size(dat.F, 2));

% Add a comment in dat
dat.dat = sprintf('Patients. Number of patients = %d (Custom colormap)', length(unique_vals));


%% Color Code by Selected Brain Region
% Read COREG annotation table
coreg = readtable('/Users/annesophiedubarry/Library/CloudStorage/SynologyDrive-NAS/0_projects/in_progress/EXPRESSONINTRA/data/COREG_test_mots_cles_amygdale_insula_post_insula_ant.tsv', ...
    "FileType", "text", 'Delimiter', '\t');

labels = 'localisation_lectrodes';  % Column containing region labels

% Define keywords to search for (only regions containing these words will be included in the color map)
region = {'insula'};  % if empty: all regions will be included

% Keep only subjects that are in the BST database
coreg(~ismember(coreg{:,1}, {chan.Channel.Comment}), :) = [];

% Find channels present both in chan and the coreg table
coreg_names = strcat(coreg.Subject, '_', coreg.Contact);
chan_names = {chan.Channel.Name};

[lia, lib] = ismember(coreg_names, chan_names);

% Check and display any unmatched electrodes
if any(lib == 0)
    warning('Some coreg electrodes were not found in chan!');
    disp(coreg_names(lib == 0));
end

% Display number of matched electrodes
fprintf('Number of matched electrodes: %d out of %d\n', sum(lia), length(coreg_names));

% Select only matching rows from coreg
coreg_filtered = coreg(lia, :);
lib_filtered = lib(lia);

% Get region labels for matched electrodes
labels_in_chan = coreg_filtered{:, labels};

% Get unique anatomical region labels and their indices (preserve original order)
[c1, ~, ic1] = unique(labels_in_chan, 'stable');

% Display number of unique regions
fprintf('There are %d unique regions found\n', length(c1));

% If no keywords defined, select all regions
if isempty(region)
    mask_region = true(size(c1));
else
    mask_region = any(cell2mat(cellfun(@(r) contains(c1, r, 'IgnoreCase', true), region, 'UniformOutput', false)), 2);
end

% Filter labels and indices based on selection
filtered_c1 = c1(mask_region);
filtered_idx = find(mask_region);

% Display selected labels
fprintf('Selected labels for color map: %s\n', strjoin(filtered_c1, ', '));

% Initialize dat.F
dat.F = zeros(length(chan.Channel), 10);

% Map region labels to numeric values (0 by default)
label_map = zeros(length(c1), 1);
label_map(filtered_idx) = 1:length(filtered_idx);  % Selected regions → 1:N

% Apply mapping to each matched electrode
values_coreg = label_map(ic1);  % size = [n_matched_electrodes x 1]

% Initialize full vector for all chan.Channel
values_chan = zeros(length(chan.Channel), 1);

% Assign values in chan space using lib_filtered
values_chan(lib_filtered) = values_coreg;

% Fill dat.F: copy region value to all time points
dat.F = repmat(values_chan, 1, size(dat.F, 2));

% Add a comment in dat
if isempty(region)
    dat.Comment = sprintf('ALL regions (%d total)', length(filtered_c1));
else
    dat.Comment = sprintf('Regions. %d regions: %s', length(filtered_c1), strjoin(filtered_c1, ', '));
end

% Set visibility of all IntraElectrodes to false
[chan.IntraElectrodes.Visible] = deal(false);

% Get names of active electrodes (those with a non-zero region value)
active_names = {chan.Channel(values_chan ~= 0).Name};
intra_names = {chan.IntraElectrodes.Name};

% Check which intra_names are matched by at least one active_name
intra_match = cellfun(@(intra) any(contains(active_names, intra)), intra_names);

% Find indices of IntraElectrodes matching at least one active electrode
idx_matched_intra = find(intra_match);

% Set visibility to true for matched IntraElectrodes
[chan.IntraElectrodes(idx_matched_intra).Visible] = deal(true);
