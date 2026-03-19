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
% Copyright (C) 2026 CNRS - Universite Aix-Marseille
% ========================================================================

function [] = mia_legacy_display_bridge(m_table_all, group_dir, list_studies_handle)
% mia_legacy_display_bridge: Open the new ROI display from a legacy study.
%
% This adapter keeps the legacy study/atlas pipeline intact:
%   legacy study MAT file -> mia_create_table_of_rois -> mia_group_gui

study_names = get(list_studies_handle, 'String');
study_idx = get(list_studies_handle, 'Value');

if isempty(study_names)
    warndlg('You must create a study first');
    return;
end

if isempty(study_idx)
    warndlg('You must select a study first');
    return;
end

if numel(study_idx) > 1
    warndlg('The new display currently opens one study at a time. Using the first selected study.');
    study_idx = study_idx(1);
end

selected_study = study_names{study_idx};
study_mat = local_get_study_file(group_dir, selected_study);

if isempty(study_mat)
    errordlg(sprintf('No study MAT file found for %s.', selected_study), 'Display Error');
    return;
end

options.signifmode = 0;
options.allow_flipsign = 0;
options.flip_thresh = -0.8;

[rois, ~, ~, ~, ~, ~, ~] = mia_create_table_of_rois(m_table_all, study_mat, options);

if isempty(rois)
    warndlg(sprintf('No ROI available to display for %s.', selected_study));
    return;
end

mia_group_gui(rois, selected_study);

end

function study_mat = local_get_study_file(group_dir, selected_study)

study_mat = '';
study_dir = fullfile(group_dir, selected_study);

if ~exist(study_dir, 'dir')
    return;
end

d = dir(fullfile(study_dir, [selected_study '*.mat']));
d = d(~[d.isdir]);

if isempty(d)
    return;
end

study_mat = fullfile(study_dir, d(1).name);

end
