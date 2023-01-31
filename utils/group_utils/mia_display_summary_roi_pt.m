function  [labels_o,colorm] = mia_display_summary_roi_pt(roi,OPTIONS) 
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
%
% 2015/5/15 : ASD creation

% Font size
FONTSZ = 12 ; 

tt=1;
for kk=1:length(roi)

    if sum(strcmp(roi{kk}.namePt,OPTIONS.ptKey))~=0
        signmoy(:,tt)= roi{kk}.signmoy(:,find(strcmp(OPTIONS.ptKey,roi{kk}.namePt)));
        labels_roi{tt} = roi{kk}.name;
        tt=tt+1;
    end
        
end

if roi{1}.freq(1)==0
    OPTIONS.title=sprintf('%s\n(LFP signal)',OPTIONS.ptKey);
else
    OPTIONS.title=sprintf('%s\nfreq band [%sH;%sHz] (step %s)',OPTIONS.ptKey,num2str(roi{1}.freq(1)),num2str(roi{1}.freq(end)),num2str(roi{1}.freq(2)-roi{1}.freq(1)));
end

% Loop over sensors
for ii=1:size(signmoy,2)

    % Finds rising and falling edges 
    tmp_onset = find(diff(abs(signmoy(:,ii))>OPTIONS.thresh)>0,1) ; 
    if isempty(tmp_onset)
        onset(ii) = size(signmoy,1);
    else
        onset(ii)= tmp_onset ;
    end
end

[Y,I] = sort(onset);

% reorder signals?labels by onsets
signmoy_o = signmoy(:,I) ; 
labels_o = labels_roi(I) ;

% set colors to timeseries 
colorm = jet(length(labels_roi)) ; 

nleft = mod(length(labels_roi),OPTIONS.nsub);
ct = 1 ;

% % sans les rois
figure('Units','normalized','position',[0 0.5 0.5 0.7],'Name',sprintf('%s Mean ROIS',OPTIONS.title),'NumberTitle','off','Position',[0.01, 0.3, 0.3, 0.3]); 
%,'Position',[0.3, 0.3, 0.3, 0.3]

maxis = ceil(max(max(abs(signmoy_o))));

for jj=1:OPTIONS.nsub

    nsig = fix(length(labels_roi)/OPTIONS.nsub);
    
    % Adds all residual timeseries on the latest graph 
    if jj==OPTIONS.nsub
        nsig = nsig +nleft ; 
    end
    % Adds residual timeseries on the first graphs (one per graph) 
%     if (nleft+1-jj)>0 
%         nsig = nsig +1 ; 
%     end
    
    set(0,'DefaultAxesColorOrder',colorm(ct:ct+nsig-1,:));
    subplot(OPTIONS.nsub,1,jj); 
    plot(roi{1}.t,signmoy_o(:,ct:ct+nsig-1)','LineWidth',1.05); legend(labels_o(ct:ct+nsig-1)) ; grid on ; ylim([-maxis maxis]); xlim([roi{1}.t(1), roi{1}.t(end)]);
%     mia_stackplot_ASD(roi{1}.t,signmoy_o(:,ct:ct+nsig-1),[1 1 1],labels_o(ct:ct+nsig-1)) ; hold on ; xlim([roi{1}.t(1), 1])%roi{1}.t(end)]);
    ct = ct+nsig ;
    title(strrep(OPTIONS.title,'_','-'),'FontSize',FONTSZ) ;
%     title(titles{jj},'FontSize',FONTSZ) ;
    xlabel('Time(s)','FontSize',FONTSZ);
    ylabel('zscore','FontSize',FONTSZ);
end

set(gcf,'color','white')
