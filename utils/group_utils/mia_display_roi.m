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

% Font size
OPTIONS.FONTSZ = 12 ; 

% Limit of Y axis (and colorbar for significance)
OPTIONS.YLIM = 1.5; 

% Default plots background color
OPTIONS.BCKGDCOLOR = [0.8275, 0.8275, 0.8275] ; 

%  Loop through regions
for ii=1:length(roi) 

    croi = roi{ii};
    
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
       hfig(ii) = figure('Name',figure_name,'Units','Pixels','Position', [1   271   432   800],'NumberTitle','off');
    else 
        % Displays one figure per ROI 
        hfig(ii) = figure('Name',figure_name,'Units','Pixels','Position',[721   271   432   800],'NumberTitle','off'); 
        % hfig = figure('Name',croi.name,'Units','Normalized','Position',[0.5 0.1 0.2  0.4],'NumberTitle','off'); 
    end
    
    % window background white 
    set(gcf,'color','white');
    
    %% Plot mean signals
    hplot = subplot(5,1,1) ;   

    for jj=1:length(croi.idPt)  
        % Plot patient averaged timeseries
        plot(croi.t,croi.signmoy(:,jj)','color',OPTIONS.clr(croi.idPt(jj),:), 'LineWidth',2); hold on ;  
        
        % % Plot halo to show variance (mean absolute deviation)
        % hPatch = plotHaloPatchMAD(hplot, croi.t, croi.signmoy(:,jj), OPTIONS.clr(croi.idPt(jj),:)) ;
        % 
    end    
    

    % Add patient' id in the legend (if they are present in the plot!)
    str_legend = arrayfun(@(i) ...
    sprintf('%s - (n = %d)', strrep(croi.namePt{i}, '_', '-'), ...
    sum(contains(croi.labels, croi.namePt{i}))), ...
    1:numel(croi.namePt), 'UniformOutput', false);

    hleg = legend(str_legend,'Location','NorthWest');


    % Add title (name of the ROI and frequencies explored (or LFP) 
    title(strrep(tmp,'_','\_'),'FontSize', OPTIONS.FONTSZ); grid on ; 
    xlim(OPTIONS.win_noedges);  ylim([-OPTIONS.YLIM ,OPTIONS.YLIM]); 
    hcol1 = colorbar ; set(hcol1,'visible','off')
    ylabel('zscore','FontSize',OPTIONS.FONTSZ);
    set(gca,'color',OPTIONS.BCKGDCOLOR);
    set(gca,'FontSize',OPTIONS.FONTSZ);
    
    % Save x-tick to apply to second planel display (masks)
    xticks =  get(gca, 'XTick') ; 
    
    %% Plot masks
    himage = subplot(5,1,[2 5]) ; 
    
    if isfield(croi,'F')
         hcol2 = display_individual_channels_image(croi.F, croi,OPTIONS) ; 
    else
        hcol2 = display_significance_masks(croi.Fmask,croi,OPTIONS) ; 
    end
  
    pos_fig = get(hfig(ii), 'position'); 
  
    %% Moves legend to the left (everything else to the rigth) 
    set(hleg,'units','pixels');
    set(hplot,'units','pixels');
    set(himage,'units','pixels');
    set(hcol1,'units','pixels');
    set(hcol2,'units','pixels');
    pos_hl = get(hleg, 'Position');
    pos_hp = get(hplot, 'position');
    pos_im = get(himage, 'position'); 

    n = 150 ;
    
    % % Write the name of the condition(s) below the legend
    % htext = text(1,1,'toto','FontSize',14) ; 
    % set(htext,'units','pixels');
    % 
    set(hfig(ii),'Position',[pos_fig(1),pos_fig(2)-n,pos_fig(3)+n,pos_fig(4)]);
    set(hplot,'position',[pos_hp(1)+n,pos_hp(2),pos_hp(3),pos_hp(4)]);
    set(himage,'position',[pos_im(1)+n,pos_im(2),pos_im(3),pos_im(4)]);
    set(hleg,'Position',[10,pos_hl(2),pos_hl(3),pos_hl(4)]);
    % set(htext,'Position',[10,pos_hl(2),0]);
    
    % Position verticale fixe en haut à droite (dans la figure)
    btnY = 960; % adapte cette valeur selon ta résolution

    % Bouton +
    btnPlus = uicontrol('Parent', hfig(ii), ...
        'Style', 'pushbutton', ...
        'String', '+', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Units', 'pixels', ...
        'Position', [hcol1.Position(1), hcol1.Position(2)+20, 20, 20], ...
        'Callback', @(src,~) adjustYLim(hplot,himage, -0.5));
    
    % Bouton -
    btnMinus = uicontrol('Parent', hfig(ii), ...
        'Style', 'pushbutton', ...
        'String', '-', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Units', 'pixels', ...
        'Position', [hcol1.Position(1), hcol1.Position(2), 20, 20], ...
        'Callback', @(src,~) adjustYLim(hplot,himage, +0.5));



end

% Set output if the function was called with an outupt variable
if nargout >= 1
    varargout{1} = hfig;
end

end

function [hcol2] = display_individual_channels_image(Fdisp, croi,OPTIONS)

    h= imagesc(croi.t,1:size(Fdisp,2),Fdisp'); caxis([-OPTIONS.YLIM OPTIONS.YLIM ]); 
    hcol2 = colorbar; xlim(OPTIONS.win_noedges)
    colormap(jet);
    
    % Light grey color for middle of the scale 
    cmap = colormap ; cmap(fix(length(cmap)/2),:) = OPTIONS.BCKGDCOLOR ; colormap(cmap) ;
    
    set(gca,'YTick',1:size(croi.Fmask,2),'YTickLabel', strrep(croi.labels,'_','\_'),'XTick',xticks); grid on ;
    % Change fontsize spearately for xtick and ytick 
    ax = gca;
    ax.XAxis.FontSize = OPTIONS.FONTSZ;
    ax.YAxis.FontSize = 10; %Smaller font for channel labels (which can be long)
    
    title(sprintf('[R_p = %0.3f]  [R_c = %0.3f]',croi.corrPt,croi.corrChan), 'FontSize', OPTIONS.FONTSZ);
    xlabel('Time(s)','FontSize',OPTIONS.FONTSZ);

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

function adjustYLim(ax, himage, delta)
    yL = ylim(ax);
    newLim = max(abs(yL)) + delta;
    newLim = max(newLim, 0.5);  % limite minimale
    ylim(ax, [-newLim newLim]);
    caxis(himage,[-newLim newLim]);
end
