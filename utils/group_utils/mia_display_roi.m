function [] = mia_display_roi(roi,OPTIONS) 
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
% Copyright (C) 2016-2021 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)
% 2015/5/15 : ASD creation

% Font size
FONTSZ = 12 ; 
BackgroundColor = [0.8275, 0.8275, 0.8275] ; 

%  Loop through regions
for ii=1:length(roi) 

    croi = roi{ii};
    
    %display frequency band
    if croi.freq(1)==0
        tmp=sprintf('%s\n(LFP signal)',croi.name);
    else
        tmp=sprintf('%s\nfreq band [%sH;%sHz] (step %sHz)',croi.name,num2str(croi.freq(1)),num2str(croi.freq(end)),num2str(croi.freq(2)-croi.freq(1)));
    end
    
    
    if strcmp(croi.name(1),'L')
        % Displays one figure per ROI 
       hfig = figure('Name',croi.name,'Units','Pixels','Position', [1   271   432   450],'NumberTitle','off');
    else 
        % Displays one figure per ROI 
        hfig = figure('Name',croi.name,'Units','Pixels','Position',[721   271   432   450],'NumberTitle','off'); 
    end
    
    % window background white 
    set(gcf,'color','white');
    
    %% Plot mean signals
    hplot = subplot(2,1,1) ;   
    
    for jj=1:length(croi.idPt)  
        % Plot patient averaged timeseries
        plot(croi.t,croi.signmoy(:,jj)','color',OPTIONS.clr(croi.idPt(jj),:), 'LineWidth',2); hold on ;  
    end    
    
    % Add patient' id in the legend (if they are present in the plot!)
    hleg = legend(strrep(croi.namePt,'_','-'),'Location','NorthWest');

    % Add title (name of the ROI and frequencies explored (or LFP) 
    title(strrep(tmp,'_','\_'),'FontSize', FONTSZ); grid on ; 
    xlim(OPTIONS.win_noedges);  ylim([-15,15]); 
    hcol1 = colorbar ; set(hcol1,'visible','off')
    ylabel('zscore','FontSize',FONTSZ);
    set(gca,'color',BackgroundColor);
    
    %% Plot masks
    himage = subplot(2,1,2) ; h= imagesc(croi.t,1:size(croi.Fmask,2),croi.Fmask'); caxis([-12 12]); hcol2 = colorbar; xlim(OPTIONS.win_noedges)
    colormap(jet);
    % Light grey color for middle of the scale 
    cmap = colormap ; cmap(33,:) = BackgroundColor ; colormap(cmap) ;
    
    set(gca,'YTick',1:size(croi.Fmask,2),'YTickLabel', strrep(croi.labels,'_','\_'),'Fontsize',8);  grid on ;
    title(sprintf('[R_p = %0.3f]  [R_c = %0.3f]',croi.corrPt,croi.corrChan), 'FontSize', FONTSZ);
    xlabel('Time(s)','FontSize',FONTSZ);

    % Moves legend to the left (everything else to the rigth) 
    set(hleg,'units','pixels');
    set(hplot,'units','pixels');
    set(himage,'units','pixels');
    set(hcol1,'units','pixels');
    set(hcol2,'units','pixels');
    pos_hl = get(hleg, 'Position');
    pos_hp = get(hplot, 'position');
    pos_im = get(himage, 'position'); 

    n = 150 ;
    pos_fig = get(hfig, 'position'); 
    
    set(hfig,'Position',[pos_fig(1),pos_fig(2)-n,pos_fig(3)+n,pos_fig(4)]);
    set(hplot,'position',[pos_hp(1)+n,pos_hp(2),pos_hp(3),pos_hp(4)]);
    set(himage,'position',[pos_im(1)+n,pos_im(2),pos_im(3),pos_im(4)]);
    set(hleg,'Position',[10,pos_hl(2),pos_hl(3),pos_hl(4)]);
    
end


