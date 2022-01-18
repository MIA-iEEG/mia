function [] = mia_display_roi_and_raster(roi,OPTIONS) 
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
        
    %% Get current roi signals 
    csig = croi.F;
    mat_csig = squeeze(cat(3,csig{:})) ; 

    % If rts were loaded
    if isfield(croi,'rts')

        % Get current roi RT 
        crt = croi.rts ; 

        % Init the reordering vector for raster plot
        I = zeros(1,size(mat_csig,2));
        ct = 1; 

     % Reorder each contact by rt
        for kk=1:size(croi.labels,1)
            % If Rts are present re-order signal 
            if ~isempty(crt{kk})

                [Yrt,Irt] = sort(cell2mat(crt{kk}),1,'descend');
                I(ct:ct+size(crt{kk},1)-1) = Irt + ct -1 ; 
                ct = ct + size(crt{kk},1) ; 
            else
                I(ct:ct+size(csig{kk},3)-1) = ct:ct+size(csig{kk},3)-1;
                ct = ct + size(csig{kk},3) ; 
            end
        end       

    else
        I = 1:size(mat_csig,2) ; 
    end

   % Displays a figure for rasters only 
    if strcmp(croi.name(end),'L')
        % Displays one figure per ROI 
       hfigrast = figure('Name',croi.name,'Units','Normalized','Position',[0 , 0.5, 0.3,1]); 
       
    else 
        % Displays one figure per ROI 
        hfigrast = figure('Name',croi.name,'Units','Normalized','Position',[0.5 , 0.5, 0.3,1]); 
    end
    imagesc(Time,1:size(mat_csig,2),mat_csig(:,I)') ;
    title(croi.name,'FontSize', FONTSZ); 
 
    c = colorbar ; caxis([-1.5,1.5]);
    c.Label.String = 'Z-Score';
    
    interticks = size(mat_csig,2)/length(croi.labels) ; 
    yticks(floor((1:length(croi.labels))*interticks)-interticks/2) ; 
    yticklabels(strrep(croi.labels,'_','\_'));
    xlabel('Time (sec)');
%     ytickangle(45) ;
    % Holds on for plotting the lines between contacts and RTs (if any)
    hold on ; 
    ax = gca;
    set(ax,'Fontsize',FONTSZ);
    % Display black dots
    ct = 0 ;
    for ii=1:size(csig,1)
        ct = size(csig{ii},3)+ ct ;
        line([Time(1), Time(end)],[ct, ct], 'color','k', 'LineWidth',2);

        % Display Rts if any  
        if isfield(croi,'rts')
            if ~isempty(crt{ii})
                current_rt = cell2mat(crt{ii})/1000 ;
                plot(sort(current_rt,'descend'),ct-length(current_rt)+1:ct,'k.','MarkerSize',12) ;
            end
        end
    end

%     colormap(jet);
%     caxis([-2,2]);
%     
    hold off ;

end

set(gcf,'color','white');
