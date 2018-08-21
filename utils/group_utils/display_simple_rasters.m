function [] = display_roi_and_raster(roi,OPTIONS) 
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
% Copyright (C) 2016-2018 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)
%
% 2015/5/15 : ASD creation

%  Loop through regions
for ii=1:length(roi) 

    croi = roi{ii};
    
    if strcmp(croi.name(end),'L')
        % Displays one figure per ROI 
       hfig = figure('Name',croi.name,'Units','Normalized','Position',[0 , 0.5, 0.3,1],'NumberTitle','off','Position',[0.01, 0.03, 0.3, 0.7]); 
    else 
        % Displays one figure per ROI 
        hfig = figure('Name',croi.name,'Units','Normalized','Position',[0.5 , 0.5, 0.3,1],'NumberTitle','off','Position',[0.01, 0.03, 0.3, 0.7]); 
        
    end
    
    % window background white 
    set(gcf,'color','white');
    
    %get back missing patient
    if ~isempty(roi{ii}.missingPt)
        tmp='';
        for ss=1:length(roi{ii}.missingPt)
            tmp=[tmp,'  ',roi{ii}.missingPt{ss}];
        end
        missingpt=sprintf('%s : missing Response Time files for %s',croi.name,tmp);
        set(gcf,'Name',missingpt);
    else
        set(gcf,'Name',croi.name);
    end 
    
    
    %% Plot mean signals
    ah1 = subplot(5,1,1) ;   
%     set(gca, 'ColorOrder',OPTIONS.clr(croi.idPt,:))
    
    for jj=1:length(croi.idPt)  
        % Plot patient averaged timeseries
        plot(croi.t,croi.signmoy(:,jj)','color',OPTIONS.clr(croi.idPt(jj),:), 'LineWidth',1.03); hold on ;  
    end    
    
    % Add patient' id in the legend (if they are present in the plot!)
    hleg = legend(strrep(croi.namePt,'_','-'),'Location','NorthWest');

    title(croi.name,'FontSize', 30); grid on ;
    grid on ; xlim(OPTIONS.win_noedges);   colorbar ; 
    ylabel('zs mean','FontSize',12);
    
    %% Plot masks
    ah2 = subplot(5,1,2) ; imagesc(croi.t,1:size(croi.Fmask,2),croi.Fmask'); caxis([-12 12]); xlim(OPTIONS.win_noedges); colorbar
    set(gca,'YTick',1:size(croi.Fmask,2),'YTickLabel', croi.labels);  grid on ;
    title(sprintf('R = %0.2f - t1 = %dms',croi.meancorr,croi.onset), 'FontSize', 28);

    %% Get current roi signals 
    csig = croi.F;
    mat_csig = squeeze(cat(3,csig{:})) ; 
        
    % Get current roi RT 
    crt = croi.rts ; 
    mat_crt = squeeze(cat(1,crt{:}));
    
    % All pt has black dots (color = 0,0,0)
    mat_colorpt = zeros(size(mat_csig,2),3);
         
    % Init the reordering vector for raster plot
    I = zeros(size(vertcat(crt{:})));
    ct = 1; 
    
    % Reorder each contact by rt
    for kk=1:size(croi.labels,1)

        [Yrt,Irt] = sort(crt{kk},1,'descend');
        I(ct:ct+size(crt{kk},1)-1) = Irt + ct -1 ; 
        ct = ct + size(crt{kk},1) ; 
    end

    %% Display the raster plot;
    
    %convert ms into s
    TIME=roi{ii}.t;
    mat_crt=(mat_crt(I)./((TIME(2)-TIME(1))*1000))*(TIME(2)-TIME(1));
         
    ah3 = subplot(5,1,[3 5]) ; 
    display_raster(roi{ii}.t, mat_csig(:,I), mat_crt, mat_colorpt(I,:),[roi{ii}.name],ah3) ;
    
    %get the right yTick and yTickLabel to display label name on raster
    %plot
        
        yTick=zeros(1,length(roi{ii}.F)*3); %code for position of yticklabel
        yTickLabel={}; %give the exact yticklabel
        sum=0;
        
        for pp=1:length(roi{ii}.F)
            
        % mean correspond to the position of label name on y axis    
        if mod(size(roi{ii}.F{pp},3),2)==0 % nb d'essai pair
            mean=size(roi{ii}.F{pp},3)./2 + sum;
        else %nb d'essai impair
            mean=(size(roi{ii}.F{pp},3)-1)./2 + sum;
        end
        
        %get back the right indexes
        ind=[1,2,3]+3*(pp-1);

        
        yTick(:,ind)=[1+sum,mean,sum+size(roi{ii}.F{pp},3)];
        yTickLabel=[yTickLabel,{'',roi{ii}.labels{pp},''}];
        
        sum=sum+size(roi{ii}.F{pp},3);
        
        end
        set(gca,'YTick',yTick)
        set(gca,'YTickLabel',yTickLabel)
 
    


    %% Adjust axes positions 
    % find window current position [x,y,width,height]
    pos3 = get(ah3,'Position');
    pos2 = get(ah2,'Position');
    pos1 = get(ah1,'Position');

    % set width of second axes equal to first
    pos3(3) = pos1(3);
    pos3(1) = pos1(1);
    set(ah3,'Position',pos3)


end

































