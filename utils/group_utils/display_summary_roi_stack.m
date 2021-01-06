function  [labels_o,colorm] = display_summary_roi_stack(roi,OPTIONS) 
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
%
% 2017/1/11 : ASD creation

% Font size
FONTSZ = 12; 

for kk=1:length(roi)

    signmoy(:,kk)= mean(roi{kk}.signmoy,2);
    labels_roi{kk} = roi{kk}.name;
end

if roi{1}.freq(1)==0
    OPTIONS.title='(LFP signal)';
else
    OPTIONS.title=sprintf('freq band [%sH;%sHz] (step %s)',num2str(roi{1}.freq(1)),num2str(roi{1}.freq(end)),num2str(roi{1}.freq(2)-roi{1}.freq(1)));
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
%figure('Units','normalized','position',[0 0.5 0.5 0.7]); 
figure('Units','normalized','position',[0 0 0.5 0.7],'Name','Mean ROIS','NumberTitle','off'); 

maxis = ceil(max(max(abs(signmoy_o))));

for jj=1:OPTIONS.nsub

    nsig = fix(length(labels_roi)/OPTIONS.nsub);
    
    if (nleft+1-jj)>0 
        nsig = nsig +1 ; 
    end
    
    set(0,'DefaultAxesColorOrder',colorm(ct:ct+nsig-1,:));
    stackplot_ASD(roi{1}.t,signmoy_o,[1 1 1],labels_o) ; hold on ; 
    
% subplot(OPTIONS.nsub,1,jj); plot(roi{1}.t,signmoy_o(:,ct:ct+nsig-1)','LineWidth',1.05); legend(labels_o(ct:ct+nsig-1)) ; grid on ; ylim([-maxis maxis]) ; xlim([roi{1}.t(1), roi{1}.t(end)]);
%     ct = ct+nsig ;
%     title(strrep(OPTIONS.title,'_','-'),'FontSize',FONTSZ) ;
%     xlabel('Time(s)','FontSize',FONTSZ);
%     ylabel('zscore','FontSize',FONTSZ);
%     
end

title(sprintf('Mean timeseries (Z-scores) (ordered at %0.2f)',OPTIONS.thresh));
    
set(gcf,'color','white')
