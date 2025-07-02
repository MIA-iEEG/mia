function varargout = mia_display_roi(roi, OPTIONS)
%
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
% Copyright (C) 2016-2022 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)
% 2015/5/15 : ASD creation
%

% Set default font size if not provided
if ~isfield(OPTIONS, 'FONTSZ')
    OPTIONS.FONTSZ = 12;
end

% Set default Y axis limits if not provided
if ~isfield(OPTIONS, 'YLIM')
    OPTIONS.YLIM = 1.5;
end

% Set default background color if not provided
if ~isfield(OPTIONS, 'BCKGDCOLOR')
    OPTIONS.BCKGDCOLOR = [0.8275, 0.8275, 0.8275];
end

% Loop through all ROI regions and plot
for ii = 1:length(roi)
    
    croi = roi{ii};
    
    % Prepare frequency band description or LFP signal label
    if croi.freq(1) == 0
        tmp = sprintf('%s\n(LFP signal)', croi.name);
    else
        stepFreq = croi.freq(2) - croi.freq(1);
        tmp = sprintf('%s\nfreq band [%sHz;%sHz] (step %sHz)', croi.name, ...
            num2str(croi.freq(1)), num2str(croi.freq(end)), num2str(stepFreq));
    end
    
    % Define figure name based on OPTIONS.Condition if present
    if isfield(OPTIONS, 'Condition')
        figure_name = OPTIONS.Condition;
    else
        figure_name = croi.name;
    end
    
    % Create figure at different positions depending on ROI name
    if strcmp(croi.name(1), 'L')
        % Left side figure position
        hfig(ii) = figure('Name', figure_name, 'Units', 'Pixels', ...
            'Position', [1 271 432 800], 'NumberTitle', 'off');
    else
        % Right side figure position
        hfig(ii) = figure('Name', figure_name, 'Units', 'Pixels', ...
            'Position', [721 271 432 800], 'NumberTitle', 'off');
    end
    
    % Set figure background color to white
    set(hfig(ii), 'Color', 'white');
    
    %% Plot mean signals in the first subplot
    hplot = subplot(5,1,1);
    hold(hplot, 'on');
    
    % Plot each patient's averaged time series with color from OPTIONS.clr
    for jj = 1:length(croi.idPt)
        plot(croi.t, croi.signmoy(:, jj)', 'Color', OPTIONS.clr(croi.idPt(jj), :), 'LineWidth', 2);
    end
    
    % Create legend with patient IDs and counts of labels
    str_legend = arrayfun(@(i) ...
        sprintf('%s - (n = %d)', strrep(croi.namePt{i}, '_', '-'), ...
        sum(contains(croi.labels, croi.namePt{i}))), ...
        1:numel(croi.namePt), 'UniformOutput', false);
    
    hleg = legend(str_legend, 'Location', 'NorthWest');
    
    % Set title, grid, axis limits, and labels
    title(strrep(tmp, '_', '\_'), 'FontSize', OPTIONS.FONTSZ);
    grid on;
    if isfield(OPTIONS, 'win_noedges')
        xlim(OPTIONS.win_noedges);
    end
    ylim([-OPTIONS.YLIM, OPTIONS.YLIM]);
    
    % Add (hidden) colorbar for alignment purposes
    hcol1 = colorbar;
    set(hcol1, 'Visible', 'off');
    
    ylabel('zscore', 'FontSize', OPTIONS.FONTSZ);
    
    % Set axes background and font size
    set(gca, 'Color', OPTIONS.BCKGDCOLOR);
    set(gca, 'FontSize', OPTIONS.FONTSZ);
    
    % Save x-tick positions for use in mask display
    xticks = get(gca, 'XTick');
    
    %% Plot significance masks or individual channels image in lower subplots
    himage = subplot(5,1, [2 5]);
    
    if isfield(croi, 'F')
        hcol2 = display_individual_channels_image(croi.F, croi, OPTIONS, xticks);
    else
        hcol2 = display_individual_channels_image(croi.Fmask, croi, OPTIONS, xticks);
    end
    
    % Get figure position for later adjustments
    pos_fig = get(hfig(ii), 'Position');
    
    % Adjust layout: Move legend to left and plots to the right
    set(hleg, 'Units', 'pixels');
    set(hplot, 'Units', 'pixels');
    set(himage, 'Units', 'pixels');
    set(hcol1, 'Units', 'pixels');
    set(hcol2, 'Units', 'pixels');
    
    pos_hl = get(hleg, 'Position');
    pos_hp = get(hplot, 'Position');
    pos_im = get(himage, 'Position');
    
    n = 150; % shift amount in pixels
    
    % Resize and reposition figure and axes
    set(hfig(ii), 'Position', [pos_fig(1), pos_fig(2)-n, pos_fig(3)+n, pos_fig(4)]);
    set(hplot, 'Position', [pos_hp(1)+n, pos_hp(2), pos_hp(3), pos_hp(4)]);
    set(himage, 'Position', [pos_im(1)+n, pos_im(2), pos_im(3), pos_im(4)]);
    set(hleg, 'Position', [10, pos_hl(2), pos_hl(3), pos_hl(4)]);
    
    % Vertical position fixed (unused variable btnY removed)
    
    % Add zoom in/out buttons for Y axis adjustment next to the colorbar
    btnPlus = uicontrol('Parent', hfig(ii), ...
        'Style', 'pushbutton', ...
        'String', '+', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Units', 'pixels', ...
        'Position', [hcol1.Position(1), hcol1.Position(2)+20, 20, 20], ...
        'Callback', @(~, ~) adjustYLim(hplot, himage, -0.5));
    
    btnMinus = uicontrol('Parent', hfig(ii), ...
        'Style', 'pushbutton', ...
        'String', '-', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Units', 'pixels', ...
        'Position', [hcol1.Position(1), hcol1.Position(2), 20, 20], ...
        'Callback', @(~, ~) adjustYLim(hplot, himage, +0.5));
    
end

% Return figure handles if requested
if nargout >= 1
    varargout{1} = hfig;
end

end


function hcol2 = display_individual_channels_image(Fdisp, croi, OPTIONS, xticks)
% Displays individual channel data as an image with colorbar and labels

imagesc(croi.t, 1:size(Fdisp, 2), Fdisp');
caxis([-OPTIONS.YLIM OPTIONS.YLIM]);
hcol2 = colorbar;
xlim(OPTIONS.win_noedges);
colormap(jet);

% Modify colormap to set middle color to background color for better visualization
cmap = colormap;
midIdx = fix(length(cmap)/2);
cmap(midIdx, :) = OPTIONS.BCKGDCOLOR;
colormap(cmap);

% Set y-axis ticks and labels (channels)
set(gca, 'YTick', 1:size(croi.Fmask, 2), 'YTickLabel', strrep(croi.labels, '_', '\_'), 'XTick', xticks);
grid on;

% Adjust font sizes for axes
ax = gca;
ax.XAxis.FontSize = OPTIONS.FONTSZ;
ax.YAxis.FontSize = 10; % Smaller font for channel labels

% Add title and x-label
title(sprintf('[R_p = %0.3f]  [R_c = %0.3f]', croi.corrPt, croi.corrChan), 'FontSize', OPTIONS.FONTSZ);
xlabel('Time(s)', 'FontSize', OPTIONS.FONTSZ);

end


function hPatch = plotHaloPatchMAD(hAxes, vTime, Y, color)
% Plot halo patch showing mean ± MAD (median absolute deviation)
%
% This function may be unused in current mia_display_roi but kept for completeness.

% Compute median absolute deviation (MAD)
madVal = median(abs(Y - repmat(median(Y', 2)', size(Y, 2), 1)')');

% Calculate mean, upper and lower bounds
avg = mean(Y, 2);
Lhi = avg + madVal';
Llow = avg - madVal';

% Draw halo patch
hPatch = patch([vTime fliplr(vTime)], [Llow' fliplr(Lhi')], [0.6 0.6 0.6], ...
    'FaceAlpha', 0.25, ...
    'FaceColor', color, ...
    'EdgeColor', 'none', ...
    'Parent', hAxes);

% Exclude patch from legend
hPatch.Annotation.LegendInformation.IconDisplayStyle = 'off';

end


function adjustYLim(ax, himage, delta)
    % Adjusts the Y axis limits of the mean signal plot and corresponding color limits
    yL = ylim(ax);
    newLim = max(abs(yL)) + delta;
    newLim = max(newLim, 0.5);  % enforce minimum limit
    ylim(ax, [-newLim, newLim]);
    caxis(himage, [-newLim, newLim]);
end
