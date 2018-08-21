function [] = display_one_roi_and_raster(roi,rasters,idpt,m_tab,OPTIONS) 
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
%     set(gca, 'ColorOrder',OPTIONS.clr(croi.idPt,:))
    
    for jj=1:length(croi.idPt)  
        % Plot patient averaged timeseries
        plot(croi.t,croi.signmoy(:,jj)','color',OPTIONS.clr(croi.idPt(jj),:), 'LineWidth',1.03); hold on ;  
    end    
    
    % Add patient' id in the legend (if they are present in the plot!)
    hleg = legend(strrep(croi.namePt,'_','-'),'Location','NorthWest');

    title(croi.name,'FontSize', FONTSZ); grid on ;
    grid on ; xlim(OPTIONS.win_noedges);   colorbar ; 
    
    %% Plot masks
    ah2 = subplot(5,1,2) ; imagesc(croi.t,1:size(croi.Fmask,2),croi.Fmask'); caxis([-12 12]); xlim(OPTIONS.win_noedges); colorbar
    set(gca,'YTick',1:size(croi.Fmask,2),'YTickLabel', strrep(croi.labels,'_','-'));  grid on ;
    title(sprintf('[R_p = %0.3f]  [R_c = %0.3f]',croi.corrPt,croi.corrChan), 'FontSize', FONTSZ);
  

    
    %% Raster 
    idx_roi = cell2mat(m_tab(:,3))==ii;
    croi = rasters(idx_roi) ; 
    mat_croi = squeeze(cat(3,croi{:})) ; 
%         crt = rts(idx_roi) ; 
%         mat_rt = squeeze(cat(1,crt{:}));
%         ccolorpt  = idpt(idx_roi)';
%         mat_colorpt = colmap(squeeze(cat(1,ccolorpt{:})),:);
        
        % All pt has black dots (color = 0,0,0)
%         mat_colorpt = zeros(size(mat_croi,2),3);
%         
%         cop = op(idx_roi) ; 
%         mat_op = squeeze(cat(1,cop{:}));
%         cidpt = idpt(idx_roi) ; 
%         mat_idpt = squeeze(cat(1,cidpt{:}));
%         
%         
%         I = zeros(size(mat_rt)); 
%         
%         % Reorder in pt order
%         [Y, Isubj] = sort(mat_idpt,1,'ascend');
%         mat_idpt = mat_idpt(Isubj);
%         mat_op = mat_op(Isubj);
%         mat_colorpt = mat_colorpt(Isubj,:);
%         mat_croi = mat_croi(:,Isubj);
%         
%         un = unique(mat_idpt) ; 
%         
%         % Reorder by RT per patient
%         for jj=1:length(un)
%         
%               idx_pt =  mat_idpt==un(jj);
%               [Y,Ipt] = sort(mat_rt(idx_pt),1,'descend');
%               I(idx_pt) = Ipt + sum(I~=0); 
%               
%               tmp = mat_rt(idx_pt) ; 
%               
%         end
% 
%         % Display the raster plot
%         ah3 = subplot(5,1,[3 5]) ; 
%         display_raster(roi{ii}.t, mat_croi(:,I), mat_rt(I),mat_colorpt(I,:),[roi{ii}.name],ah3) ;
%         
%         

       %# find window current position [x,y,width,height]
        pos3 = get(ah3,'Position');
        pos2 = get(ah2,'Position');
        pos1 = get(ah1,'Position');

        %# set width of second axes equal to first
        pos3(3) = pos1(3);
        pos3(1) = pos1(1);
        set(ah3,'Position',pos3)


end


