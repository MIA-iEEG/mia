function [id, labels, idx] = mia_filter_localized_channels(m_table_as, idx_subjloc , id_contact, id_ncontact, bilabels, inclusive)
% ***********************************************************************
% ***********************************************************************
% Copyright (C) 2016-2022 CNRS - Universite Aix-Marseille
%
% This software was developed by
% Anne-Sophie Dubarry (CNRS University of Aix-Marseille)
%
% ***********************************************************************
%
% m_table_as : table of localizations (all patients)
% idx_subjloc : Indices of the current patient in m_table_as
% id_ncontact : Index of the contacts in m_table_as
% Bilabels : Bipolar channels present in the data
% TODO: add an OPTION for Bilpolar Channels: 
% conservative labeling (only channels where the 2 contacts are in the different ROIs are labeled with both ROI of contacts1 and contact2)
% ROI labelled ROI
% Restrictive labeling (bipolar channels are labeled with both
% ROI of contacts1 and contact2) 
% Contributors : Dewmith Weerasena

% Initialize labels as a cell array to store channel names.
% Previous versions used labels = [], which caused concatenation errors
% when adding string labels (e.g., 'A_1', 'A_1_B_2') with CAT.
labels = {};
idx = [];
id = [];

% Compatibility: allow char arrays
if ischar(bilabels)
    bilabels = cellstr(bilabels);
end

% ----------------------------------------------------------------------
% Build canonical contact labels from the localization table
% ----------------------------------------------------------------------
% The localization table stores electrode shaft and contact number separately.
% Example:
%   shaft = 'A'
%   contact = 1
% This block concatenates them into a normalized contact label: A1
contacts = strrep(strcat(m_table_as(idx_subjloc,id_contact), ...
    num2str(cell2mat(m_table_as(idx_subjloc,id_ncontact)))), ' ', '');

% ----------------------------------------------------------------------
% Loop through every channel label present in the data
% ----------------------------------------------------------------------
for jj = 1:length(bilabels)
    current_label = bilabels{jj};
    if isempty(current_label)
        continue;
    end

    % ------------------------------------------------------------------
    % MONOPOLAR MATCHING
    % ------------------------------------------------------------------
    % Some datasets store monopolar labels with underscores (e.g., B_1).
    % The localization table however uses compact labels (B1).
    %
    % Therefore we try two candidates:
    %   - the original label
    %   - the label with underscores removed
    %
    % If either matches the localization table, the channel is accepted
    % as a monopolar contact.
    mono_candidates = unique({
        current_label
        strrep(current_label, '_', '') });
    
    idx1 = find(ismember(contacts, mono_candidates));

    if ~isempty(idx1)
        n1 = numel(idx1);
        labels = cat(1, labels, repmat({current_label}, n1, 1));
        idx = cat(1, idx, idx1(:));
        id = cat(1, id, repmat(jj, n1, 1));
        continue;
    end

    % ------------------------------------------------------------------
    % BIPOLAR MATCHING
    % ------------------------------------------------------------------
    % If the label was not matched as monopolar, attempt to interpret it
    % as a bipolar channel (two contacts separated by "_").
    %
    % Supported formats include:
    %   A1_B2
    %   A_1_B_2
    %   B1_A_2
    %   B_1_A_4
    %   A'_1_B'_2
    %
    % The regex extracts the first and second contact tokens.
    tokens = regexp(current_label, '^([A-Za-z'']+(?:_?p)?_?\d+)_([A-Za-z'']+(?:_?p)?_?\d+)$', 'tokens', 'once');


    % If the label cannot be parsed as bipolar, skip it
    if isempty(tokens)
        continue;
    end

    % Extract the two contacts
    first = tokens{1};
    second = tokens{2};

    % ------------------------------------------------------------------
    % Normalize contact format
    % ------------------------------------------------------------------
    % The localization table stores contacts without underscores (A1),
    % so we remove underscores from the parsed tokens before matching.
    first_norm = strrep(first, '_', '');
    second_norm = strrep(second, '_', '');

    % Find matches for each contact
    idx1 = find(ismember(contacts, first_norm));
    idx2 = find(ismember(contacts, second_norm));

    % ------------------------------------------------------------------
    % CASE 1: both contacts are localized
    % ------------------------------------------------------------------
    % In this case the bipolar channel is duplicated so that each contact
    % is linked to its corresponding localization row (legacy behavior).
    if ~isempty(idx1) && ~isempty(idx2)
        n12 = numel(idx1) + numel(idx2);
        labels = cat(1, labels, repmat({current_label}, n12, 1));
        idx = cat(1, idx, idx1(:));
        idx = cat(1, idx, idx2(:));
        id = cat(1, id, repmat(jj, n12, 1));

    % ------------------------------------------------------------------
    % CASE 2: only first contact localized (inclusive mode)
    % ------------------------------------------------------------------
    elseif inclusive ~= 0 && ~isempty(idx1) && isempty(idx2)
        n1 = numel(idx1);
        labels = cat(1, labels, repmat({current_label}, n1, 1));
        idx = cat(1, idx, idx1(:));
        id = cat(1, id, repmat(jj, n1, 1));

    % ------------------------------------------------------------------
    % CASE 3: only second contact localized (inclusive mode)
    % ------------------------------------------------------------------
    elseif inclusive ~= 0 && isempty(idx1) && ~isempty(idx2)
        n2 = numel(idx2);
        labels = cat(1, labels, repmat({current_label}, n2, 1));
        idx = cat(1, idx, idx2(:));
        id = cat(1, id, repmat(jj, n2, 1));
    end
end

end
