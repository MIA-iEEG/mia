function [] = display_roi_and_raster(roi,rasters,idpt,m_tab,OPTIONS) 
%
% ***********************************************************************
%
%  Copyright (C) 2015 CNRS - Universite Aix-Marseille
% 
%  This software was developed by
%       Anne-Sophie Dubarry (1) (2) 
%
%       (1) CNRS Universite Aix-Marseille
%       (2) INSERM 
%
% ***********************************************************************
%
% 2015/5/15 : ASD creation

FONTSZ = 18 ;

%  Loop through regions
for ii=1:length(roi) 

    croi = roi{ii};
    Time = croi.t ; 
    
    if strcmp(croi.name(end),'L')
        % Displays one figure per ROI 
       hfig = figure('Name',croi.name,'Units','Normalized','Position',[0 , 0.5, 0.3,0.8]); 
    else 
        % Displays one figure per ROI 
        hfig = figure('Name',croi.name,'Units','Normalized','Position',[0.5 , 0.5, 0.3,0.8]); 
    end
    
    % window background white 
    set(gcf,'color','white');
    
    %% Plot mean signals
    ah1 = subplot(5,1,1) ;   
    
    for jj=1:length(croi.idPt)  
        % Plot patient averaged timeseries
        plot(Time,croi.signmoy(:,jj)','color',OPTIONS.clr(croi.idPt(jj),:), 'LineWidth',1.03); hold on ;  
    end    
    
    % Add patient' id in the legend (if they are present in the plot!)
    hleg = legend(strrep(croi.namePt,'_','-'),'Location','NorthWest');

    title(croi.name,'FontSize', FONTSZ); grid on ;
    grid on ; xlim(OPTIONS.win_noedges);   colorbar('Visible','off') ; 
    
    m = quantile2(abs(croi.signmoy(:)),.99) ; 
    ylim([-abs(m),abs(m)]);
    %% Plot masks
    ah2 = subplot(5,1,2) ; imagesc(Time,1:size(croi.Fmask,2),croi.Fmask'); caxis([-12 12]); xlim(OPTIONS.win_noedges); colorbar
    set(gca,'YTick',1:size(croi.Fmask,2),'YTickLabel', strrep(croi.labels,'_','-'));  grid on ;
    title(sprintf('[R_p = %0.3f]  [R_c = %0.3f]',croi.corrPt,croi.corrChan), 'FontSize', FONTSZ);
  
    %% Raster
    idx_roi = cell2mat(m_tab(:,3))==ii;
    croi = rasters(idx_roi) ; 
    mat_croi = squeeze(cat(3,croi{:})) ; 

    % Displays raster plots
    ah3 = subplot(5,1,[3 5]) ; 
   
    % find window current position [x,y,width,height]
    pos3 = get(ah3,'Position');
    pos2 = get(ah2,'Position');
    pos1 = get(ah1,'Position');

    % set width of second axes equal to first
    pos3(3) = pos1(3);
    pos3(1) = pos1(1);
    set(ah3,'Position',pos3)


    imagesc(ah3,Time,1:size(mat_croi,2),mat_croi') ; 
    hold on ; 

    % Display lines to separate channels
    ct = 0 ;
    for ii=1:size(croi,1)
        ct = size(croi{ii},3)+ ct ;
        plot(ah3,Time,ct,'k.') ; 

    end

    hold off ; colorbar ; caxis([-2,2]);

  
end


