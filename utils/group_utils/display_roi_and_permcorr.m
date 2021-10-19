function [] = display_roi_and_permcorr(roi,OPTIONS) 
%
% ***********************************************************************
%
%  Copyright (C) 2016-2018 CNRS - Universite Aix-Marseille
%
%  This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)
%
% ***********************************************************************
% 2015/5/15 : ASD creation

% Font size
FONTSZ = 12 ; 
BackgroundColor = [ 0.8275, 0.8275, 0.8275] ; 

fprintf('-------------------------------------------------------------------\n');
fprintf('Table description :   \n');
fprintf('Number of patients - Region : actual correlation < Quantile+5%%  \n');
fprintf('-------------------------------------------------------------------\n');
        
%  Loop through regions
for ii=1:length(roi) 

    croi = roi{ii};
    
    % If only one channel do nothing (no correlations to display)
    if length(croi.labels) <=1 ; continue ; end 
    
    %display frequency band
    if croi.freq(1)==0
        tmp=sprintf('%s\n(LFP signal)',croi.name);
    else
        tmp=sprintf('%s\nfreq band [%sH;%sHz] (step %s)',croi.name,num2str(croi.freq(1)),num2str(croi.freq(end)),num2str(croi.freq(2)-croi.freq(1)));
    end
    
    
    if strcmp(croi.name(1),'L')
        % Displays one figure per ROI 
       hfig = figure('Name',croi.name,'Units','Pixels','Position', [1   271   432   450],'NumberTitle','off');
%        hfig = figure('Name',croi.name,'Units','Normalized','Position',[0 , 0.3, 0.3,0.5],'NumberTitle','off'); 
    else 
        % Displays one figure per ROI 
        hfig = figure('Name',croi.name,'Units','Pixels','Position',[721   271   432   450],'NumberTitle','off'); 
%         hfig = figure('Name',croi.name,'Units','Normalized','Position',[0.5 , 0.3, 0.3,0.5],'NumberTitle','off'); 
    end
    
    % window background white 
    set(gcf,'color','white');
    
    %% Plot mean signals
    hplot = subplot(3,1,1) ;   
%     set(gca, 'ColorOrder',OPTIONS.clr(croi.idPt,:))
    
    for jj=1:length(croi.idPt)  
        % Plot patient averaged timeseries
        plot(croi.t,croi.signmoy(:,jj)','color',OPTIONS.clr(croi.idPt(jj),:), 'LineWidth',2); hold on ;  
    end    
    
    % Add patient' id in the legend (if they are present in the plot!)
    hleg = legend(strrep(croi.namePt,'_','-'),'Location','NorthWest');

    title(tmp,'FontSize', FONTSZ); grid on ; colorbar ;
    xlim(OPTIONS.win_noedges);   
    % ASD 2017/11/2
    ylim([-15,15]); 
    hcol1 = colorbar ; 
    ylabel('zscore','FontSize',FONTSZ);
    set(gca,'color',BackgroundColor);
    
    %% Plot masks
    himage = subplot(3,1,2) ; h= imagesc(croi.t,1:size(croi.Fmask,2),croi.Fmask'); caxis([-12 12]); hcol2 = colorbar; xlim(OPTIONS.win_noedges)
    colormap(jet);
    % Light grey color for middle of the scale 
    cmap = colormap ; cmap(33,:) = BackgroundColor ; colormap(cmap) ;
    set(gca,'YTick',1:size(croi.Fmask,2),'YTickLabel', strrep(croi.labels,'_','\_'));  grid on ;
    title(sprintf('[Rp = %0.2f] [Rc = %0.2f] [1st Onset = %0.3fsec]',croi.corrPt,croi.corrChan,croi.onset), 'FontSize', FONTSZ);
%     xlabel('Time(s)','FontSize',FONTSZ);
    
    %% Plot distrib correlation
    hhist = subplot(3,1,3) ; hist(croi.corr_permut,100); xlabel('Surrogate distribution of max correlation (permutations)');
    
%     title(sprintf('Q+5%s : %1.3f ,Q-5%s : %1.3f ','%',mia_quantile(croi.corr_permut(croi.corr_permut>0),0.95),'%',mia_quantile(croi.corr_permut(croi.corr_permut<0),0.95)));
%   
    %For NOW only look at positive values
    title(sprintf('Q+5%s : %1.3f ','%',mia_quantile(croi.corr_permut(croi.corr_permut>0),0.95)));
    
    if croi.corrChan>=mia_quantile(croi.corr_permut(croi.corr_permut>0),0.95) 
        fprintf(sprintf('%0.0f - %s : %1.3f > %1.3f *** \n',length(croi.idPt),croi.name,croi.corrChan,mia_quantile(croi.corr_permut(croi.corr_permut>0),0.95)));
    else
        fprintf(sprintf('%0.0f - %s : %1.3f < %1.3f\n',length(croi.idPt),croi.name,croi.corrChan,mia_quantile(croi.corr_permut(croi.corr_permut>0),0.95)));
    end
    % Moves legend to the left (everything else to the rigth) 
    set(hleg,'units','pixels');
    set(hplot,'units','pixels');
    set(himage,'units','pixels');
    set(hhist,'units','pixels');
    set(hcol1,'units','pixels');
    set(hcol2,'units','pixels');
    pos_hl = get(hleg, 'Position');
    pos_hp = get(hplot, 'position');
    pos_im = get(himage, 'position'); 
    pos_hi = get(hhist, 'position'); 

    n = 150 ;
    pos_fig = get(hfig, 'position'); 
    
    set(hfig,'Position',[pos_fig(1),pos_fig(2)-n,pos_fig(3)+n,pos_fig(4)]);
    set(hplot,'position',[pos_hp(1)+n,pos_hp(2),pos_hp(3),pos_hp(4)]);
    set(himage,'position',[pos_im(1)+n,pos_im(2),pos_im(3),pos_im(4)]);
    set(hhist,'position',[pos_hi(1)+n,pos_hi(2),pos_hi(3),pos_hi(4)]);
    set(hleg,'Position',[10,pos_hl(2),pos_hl(3),pos_hl(4)]);
       
 end


