function [hfig] = mia_display_roi_conditions(rois,OPTIONS) 
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

% Font size
OPTIONS.FONTSZ = 12 ; 

% Limit of Y axis (and colorbar for significance)
OPTIONS.YLIM = 1.5; 

% Default plots background color
OPTIONS.BCKGDCOLOR = [0.8275, 0.8275, 0.8275] ; 

Conditions = strsplit(OPTIONS.Condition,'-');

%  Loop through regions
for ii=1:length(rois) 

    croi = rois{ii};
    
    %display frequency band
    if croi.freq(1)==0
        tmp=sprintf('%s\n(LFP signal)',croi.name);
    else
        tmp=sprintf('%s\nfreq band [%sH;%sHz] (step %sHz)',croi.name,num2str(croi.freq(1)),num2str(croi.freq(end)),num2str(croi.freq(2)-croi.freq(1)));
    end
    
    % Define figure name
    if isfield(OPTIONS,'Condition') ; figure_name = OPTIONS.Condition ; else figure_name =croi.name; end

    if strcmp(croi.name(1),'L')
        % Displays one figure per ROI 
       hfig(ii) = figure('Name',figure_name,'Units','Pixels','Position', [1   271   600   500],'NumberTitle','off');
    else 
        % Displays one figure per ROI 
        hfig(ii) = figure('Name',figure_name,'Units','Pixels','Position',[721   271   600 500],'NumberTitle','off'); 
        % hfig = figure('Name',croi.name,'Units','Normalized','Position',[0.5 0.1 0.2  0.4],'NumberTitle','off'); 
    end
    
    % window background white 
    set(gcf,'color','white');

    
    %% Plot mean signals
    hplot = subplot(2,1,1) ;   
 
    for iCond=1:length(Conditions)  
    
        % Plot patient averaged timeseries
        plot(croi.t,mean(croi.signmoyAll{iCond},2)','color',OPTIONS.clr(iCond,:), 'LineWidth',2); hold on ;  
    end
    
    hleg2 = legend(Conditions,'Location','NorthWest');
    
    xlim(OPTIONS.win_noedges);  ylim([-OPTIONS.YLIM ,OPTIONS.YLIM]); 
    hcol1 = colorbar ; set(hcol1,'visible','off')
    ylabel('zscore','FontSize',OPTIONS.FONTSZ);
    set(gca,'color',OPTIONS.BCKGDCOLOR);
    set(gca,'FontSize',OPTIONS.FONTSZ);
    
    hcol2 = colorbar ; set(hcol2,'visible','off')
    grid on ; 
 
    %% Plot the grand averages
    hplot2 = subplot(2,1,2); 
     
    cond = [] ; 

    for iCond=1:length(Conditions)  
      
        % Plot patient averaged timeseries
        plot(croi.t,croi.signmoyAll{iCond}','color',OPTIONS.clr(iCond,:), 'LineWidth',2); hold on ;  
        
        % % Plot halo to show variance (mean absolute deviation)
        % hPatch = plotHaloPatchMAD(hplot, croi.t, croi.signmoyAll{iCond}, OPTIONS.clr(iCond,:)) ;
        
        % Assuming croi.namePt is a 4x1 cell array and Conditions is a cell array of the same length
        combined = cellfun(@(name) strcat(name, ' - ', Conditions{iCond}), croi.namePt, 'UniformOutput', false);
     
        % Optionally, concatenate the combined condition name for further processing
        cond = [cond; combined];  % Append the result to 'cond', if needed
    
    end    

    hleg = legend(cond,'Location','NorthWest');

    % Add title (name of the ROI and frequencies explored (or LFP) 
    title(strrep(tmp,'_','\_'),'FontSize', OPTIONS.FONTSZ); grid on ; 
    xlim(OPTIONS.win_noedges);  ylim([-OPTIONS.YLIM ,OPTIONS.YLIM]); 
    hcol1 = colorbar ; set(hcol1,'visible','off')
    ylabel('zscore','FontSize',OPTIONS.FONTSZ);
    set(gca,'color',OPTIONS.BCKGDCOLOR);
    set(gca,'FontSize',OPTIONS.FONTSZ);
    
    % Save x-tick to apply to second planel display (masks)
    xticks =  get(gca, 'XTick') ; 
     pos_fig = get(hfig(ii), 'position'); 
 

    %% Moves legend to the left (everything else to the rigth) 
    set(hleg,'units','pixels');
    set(hleg2,'units','pixels');
    set(hplot,'units','pixels');
    set(hplot2,'units','pixels');
    set(hcol1,'units','pixels');
    set(hcol2,'units','pixels');
    pos_hl = get(hleg, 'Position');
    pos_hp = get(hplot, 'position');
    pos_hp2 = get(hplot2, 'position'); 
    pos_hl2 = get(hleg2, 'Position');
   
    n = 150 ;
    
    % % Write the name of the condition(s) below the legend
    % htext = text(1,1,'toto','FontSize',14) ; 
    % set(htext,'units','pixels');
    % 
    set(hfig(ii),'Position',[pos_fig(1),pos_fig(2)-n,pos_fig(3)+n,pos_fig(4)]);
    set(hplot,'position',[pos_hp(1)+n,pos_hp(2),pos_hp(3),pos_hp(4)]);
    set(hplot2,'position',[pos_hp2(1)+n,pos_hp2(2),pos_hp2(3),pos_hp2(4)]);
    set(hleg,'Position',[10,pos_hl(2),pos_hl(3),pos_hl(4)]);
    set(hleg2,'Position',[10,pos_hl2(2),pos_hl2(3),pos_hl2(4)]);
    
    % Bouton +
    btnPlus = uicontrol('Parent', hfig(ii), ...
        'Style', 'pushbutton', ...
        'String', '+', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Units', 'pixels', ...
        'Position', [hcol1.Position(1), hcol1.Position(2)+20, 20, 20], ...
        'Callback', @(src,~) adjustYLim(hplot2, -0.5));
    
    % Bouton -
    btnMinus = uicontrol('Parent', hfig(ii), ...
        'Style', 'pushbutton', ...
        'String', '-', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Units', 'pixels', ...
        'Position', [hcol1.Position(1), hcol1.Position(2), 20, 20], ...
        'Callback', @(src,~) adjustYLim(hplot2, +0.5));

    % Bouton +
    btnPlus = uicontrol('Parent', hfig(ii), ...
        'Style', 'pushbutton', ...
        'String', '+', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Units', 'pixels', ...
        'Position', [hcol2.Position(1), hcol2.Position(2)+20, 20, 20], ...
        'Callback', @(src,~) adjustYLim(hplot, -0.5));
    
    % Bouton -
    btnMinus = uicontrol('Parent', hfig(ii), ...
        'Style', 'pushbutton', ...
        'String', '-', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Units', 'pixels', ...
        'Position', [hcol2.Position(1), hcol2.Position(2), 20, 20], ...
        'Callback', @(src,~) adjustYLim(hplot, +0.5));
    
end
end
function hPatch = plotHaloPatchMAD(hAxes, vTime, Y, color)

% Compute median absolute deviation (MAD)
% Note: mad = median(abs(X - median(X)))
mad = median(abs(Y-repmat(median(Y'),size(Y,2),1)')');

% Set mean, low and high halo borders
avg=mean(Y,2) ;
Lhi=mean(Y,2) + mad' ; 
Llow=mean(Y,2) - mad' ; 
    
% Draw halo patch 
hPatch = patch([vTime fliplr(vTime)], [Llow'  fliplr(Lhi')], [0.6 0.6 0.6],...
        'FaceAlpha', 0.25, ...
        'FaceColor', color, ...
        'EdgeColor', 'none', ...
        'Parent',    hAxes);

% Skip the name of the previous plot from the legend
hPatch.Annotation.LegendInformation.IconDisplayStyle = 'off';

end

function adjustYLim(ax, delta)
    yL = ylim(ax);
    newLim = max(abs(yL)) + delta;
    newLim = max(newLim, 0.5);  % limite minimale
    ylim(ax, [-newLim newLim]);
    
end
