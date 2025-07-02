function [roi, varargout] = mia_get_roi_bst(m_table_effect, t, s, smask, freqs, opt)
% mia_get_roi_bst: Build region-of-interest (ROI) data structure from Brainstorm signals.
%
% Copyright (C) 2025 Anne-Sophie Dubarry
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%
% USAGE:
%   [roi, roi_names] = mia_get_roi_bst(m_table_effect, t, s, smask, freqs, opt)
%
% INPUTS:
%   - m_table_effect : Cell array [n x 6], with subject, electrode info, hemisphere, ROI, and labels
%   - t              : Time vector
%   - s              : Data matrix [channels x time x freq]
%   - smask          : Significance mask, same size as s
%   - freqs          : Cell array of frequency labels
%   - opt            : Struct with fields:
%                      * freq        : Frequency index
%                      * nPt         : Minimum number of patients per ROI
%                      * signifmode  : 0 (all), 1 (only significant)
%                      * signmode    : 'signed' or 'abs'
%                      * montage     : 'bipolar' or 'monopolar'
%                      * flip_thresh : (optional) threshold for polarity flipping
%
% OUTPUTS:
%   - roi        : Struct array with one element per ROI
%   - roi_names  : (optional) Cell array of ROI names

% Column indices in m_table_effect
id_subj = 1;
id_loc  = 5;
id_lat  = 4;
id_label = 6;

% Prepare ROI identifiers
roi_names = strcat(m_table_effect(:,id_lat), '.', m_table_effect(:,id_loc));
[unique_rois, ~, ~] = unique(roi_names, 'stable');
[unique_subj, ~, ~] = unique(m_table_effect(:,id_subj));
nSubj = numel(unique_subj);

% Colors for visualization (not used directly here)
% clr = hsv(nSubj); % optionally retained for GUI

% Determine which contacts to include (all or only significant)
if opt.signifmode ~= 0
    signif_mask = cell2mat(m_table_effect(:, id_label + opt.freq));
else
    signif_mask = ones(size(m_table_effect,1),1);
end

roi = {};
ctRoi = 1;

% Loop over each unique ROI
for jj = 1:length(unique_rois)
    
    is_roi = ismember(roi_names, unique_rois(jj));
    in_roi = find(is_roi & signif_mask);

    if isempty(in_roi)
        continue;
    end

    % Subjects contributing to this ROI
    subj_active = unique(m_table_effect(in_roi, id_subj));
    subj_indices = find(ismember(unique_subj, subj_active));

    % Skip if not enough subjects
    if numel(subj_active) < opt.nPt
        continue;
    end

    % Keep only first occurrence per channel (remove duplicates)
    [~, ia, ~] = unique(m_table_effect(in_roi, id_label));
    idx_signals = in_roi(ia);

    % Extract data
    all_sig    = s(idx_signals, :, opt.freq)';
    masked_sig = smask(idx_signals, :, opt.freq)';
    labels_roi = m_table_effect(idx_signals, id_label);

    % Extract patient names
    splitLabels = cellfun(@(x) strsplit(x, '_'), labels_roi, 'UniformOutput', false);
    ptname = cellfun(@(x) x{1}, splitLabels, 'UniformOutput', false);

    % Flip signals if option is used
    if isfield(opt, 'flip_thresh')
        [all_sig, masked_sig, labels_roi, r] = flip_signals(all_sig, masked_sig, labels_roi, opt.flip_thresh);
    else
        r = corrcoef(all_sig);
    end

    % Compute mean signal per patient
    [~, ~, IC] = unique(ptname);
    mean_sig_subj = zeros(size(all_sig,1), max(IC));
    for ss = 1:max(IC)
        if strcmp(opt.signmode, 'signed')
            mean_sig_subj(:, ss) = mean(all_sig(:, IC == ss), 2);
        else
            mean_sig_subj(:, ss) = mean(abs(all_sig(:, IC == ss)), 2);
        end
    end

    % Compute inter-subject correlation
    rPt = corrcoef(mean_sig_subj);
    trPt = tril(atanh(rPt), -1);  % Fisher transform and extract lower triangle
    roi{ctRoi}.corrPt = tanh(mean(trPt(trPt ~= 0)));
    roi{ctRoi}.allCorrPt = trPt(trPt ~= 0);

    % Compute inter/intra-channel correlation
    trChan = tril(r, -1);
    roi{ctRoi}.corrChan = tanh(mean(trChan(trChan ~= 0)));
    roi{ctRoi}.allCorrChan = trChan(trChan ~= 0);

    % Find onset of the first significant signal
    activity = sum(masked_sig', 1);
    change_points = find(diff(activity) ~= 0, 1);
    if ~isempty(change_points)
        roi{ctRoi}.onset = t(change_points);
    else
        roi{ctRoi}.onset = NaN;
    end

    % Populate ROI fields
    roi{ctRoi}.name    = unique_rois{jj};
    roi{ctRoi}.signmoy = mean_sig_subj;
    roi{ctRoi}.Fmask   = masked_sig;
    roi{ctRoi}.labels  = labels_roi;
    roi{ctRoi}.idPt    = subj_indices;
    roi{ctRoi}.namePt  = unique_subj(subj_indices);
    roi{ctRoi}.t       = t;
    roi{ctRoi}.freq    = freqs{opt.freq};
    roi{ctRoi}.F       = all_sig;

    ctRoi = ctRoi + 1;
end

% Optional output: list of ROI names
if nargout > 1
    varargout{1} = unique_rois;
end

end

function [all_sig, masked_sig, labels_roi, r] = flip_signals(all_sig, masked_sig, labels_roi, flip_thresh)
% Flip signals if they are anti-correlated beyond threshold

    r = corrcoef(all_sig);
    A = (r < flip_thresh);
    Ap = (r > -flip_thresh);
    
    % Identify row with maximum anti-correlation
    [Y, I] = max(sum(A));
    if Y == 0
        I = [];
    end
    
    % Flip mask: +1 for keep, -1 for flip
    fl = ones(1, size(A, 1));
    if ~isempty(I)
        fl(Ap(I, :)) = -1;
    end

    % Apply flipping
    all_sig = all_sig .* fl;
    masked_sig = masked_sig .* fl;

    % Update labels
    labels_roi(fl == -1) = strcat(labels_roi(fl == -1), '_FLP');

    % Recompute correlation matrix after flipping
    r = corrcoef(all_sig);
end
