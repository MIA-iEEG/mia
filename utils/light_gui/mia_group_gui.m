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
%
% This function creates a user interface to visualize sEEG data per regions
% of interest (ROI)

function [] = mia_group_gui(varargin)
 
rois_cond = vertcat(varargin{1:nargin-1}) ; 
Condition = varargin{end}; 

% Create UIFigure
UIFigure = uifigure('Name', sprintf('MIA GAnalysis : %s',Condition), ...
                    'Units','normalized', ...
                    'Position', [0.1 0 0.29 0.9], ...
                    'HandleVisibility', 'on', ...
                    'Color', [0.9 1.0 0.9]);%, ...
                    % 'WindowStyle', 'alwaysontop');
% Set a layout manager on the whole figure
UIFigureLayout = uigridlayout(UIFigure, [1, 1]);
UIFigureLayout.RowHeight = {'1x'};
UIFigureLayout.ColumnWidth = {'1x'};

% Place the TabGroup inside the layout
TabGroup = uitabgroup(UIFigureLayout);
TabGroup.Layout.Row = 1;
TabGroup.Layout.Column = 1;

Condition = strsplit(Condition,'-');

% Prepare one tab per condition extracted
for iCond=1:size(rois_cond,1)
    
    % ---- Create one Tab per condition ----
    UITab = uitab(TabGroup, 'Title', Condition{iCond});

    % Create a main grid layout in the tab: 2 rows (table + buttons)
    MainLayout = uigridlayout(UITab, [2, 1]);
    MainLayout.RowHeight = {'1x', 50};  % Table fills, button row is 50px
    MainLayout.ColumnWidth = {'1x'};

    mia_display_table_rois(MainLayout, rois_cond(iCond,:), Condition{iCond})
    
end

% ---- Create a Tab with all group of tabs ----
UITab = uitab(TabGroup, 'Title', 'Group');

% Create a main grid layout in the tab: 2 rows (table + buttons)
MainLayout = uigridlayout(UITab, [2, 1]);
MainLayout.RowHeight = {'1x', 50};  % Table fills, button row is 50px
MainLayout.ColumnWidth = {'1x'};

mia_display_table_rois(MainLayout, rois_cond, Condition)
