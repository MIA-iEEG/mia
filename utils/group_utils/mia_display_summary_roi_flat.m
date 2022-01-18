function  [labels_o,colorm] = mia_display_summary_roi_flat(roi,OPTIONS) 
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
% 2016/10/13 : ASD creation
% Display black and white averaged time series per regions

FONTSZ = 12 ;

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
% Write out number of electrodes/patients per regions
fprintf('Number_electrode / Number of patients : \n');
for ii=1:length(roi) ; fprintf(sprintf('%s : %d/%d\n',roi{I(ii)}.name,size(roi{I(ii)}.Fmask,2),size(roi{I(ii)}.signmoy,2))) ; end

%Here this is specific to flat B&W display

figure ; imagesc(roi{1}.t,1:length(labels_o),signmoy_o'.*(abs(signmoy_o')>OPTIONS.threshdisp)) ; set(gca,'YTick',1:length(labels_o),'YTickLabel',labels_o) ; 
% toto = colormap ; toto(1,:) = [1 1 1]; colormap(toto) ; grid on ; set(gcf,'color','white');
colormap jet; %toto = flipdim(colormap,1); toto(1,:) = [1 1 1];
grid on ; set(gcf,'color','white');colorbar ;set(gca,'Fontsize',FONTSZ); xlabel('Time (sec'); title(sprintf('Mean timeseries (ordered at %0.2f zscore) ; Threshold at %0.2f',OPTIONS.thresh,OPTIONS.threshdisp));

% toto = colormap ; 
% toto(1,:) = [1 1 1];
% colormap(toto) ; colorbar ; caxis([1 7.69]) ; set(gca,'Fontsize',20);
% colormap(toto) ; colorbar ;  set(gca,'Fontsize',20);
caxis([-12 12]) ;
% 
% colorbounds = roi{1}.t(roi{1}.t>=0) ; 
% colorm = jet(size(roi{1}.t,2)) ; 
% % colorm = colorm(Y-sum(roi{1}.t<0),:) ; 
% colorm = colorm(Y,:) ; 
