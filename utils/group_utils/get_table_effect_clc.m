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
function [m_table_effect, s, smask, all_labels] = get_table_effect_clc(m_table_as, ganalysis) 

id_contact = 2; 
id_ncontact = 3;
id_subj = 1 ; 

threshp = [] ; 
m_table_effect = [] ;
idx_ganalysisloc = [] ;
all_labels = [] ;
s = [];
smask = [];
Fs = 1/(ganalysis{1}.Time(2)-ganalysis{1}.Time(1)) ; 

% Adds 'p' for left contact in labels %% for CLC
% m_table_as(strcmp(m_table_as(:,4),'L'),2) = strcat(m_table_as(strcmp(m_table_as(:,4),'L'),2),'p') ; 

% In case there is already a 'p' at the end of the electrode name (if added
% there will be no more correspondance with the data
% left_labels = m_table_as(strcmp(m_table_as(:,4),'L'),2) ;

% ASD for now loop throught the localization table ; COULD BE vectorized? 
for ii=1:size(m_table_as,1)
    % If contact is left lateralized
   if strcmp(m_table_as(ii,4),'L')
        tmp = cell2mat(m_table_as(ii,2));    
        last_char = tmp(end); 
    if ~strcmp(last_char,'p')
        % Add a 'p' at the end of the contact name if it is left
        % lateralized (CCF)
        
        m_table_as(ii,2) = strcat(m_table_as(ii,2),'p');
%         m_table_as( strcmp(m_table_as(:,4),'L'),2) = strcat(m_table_as(strcmp(m_table_as(:,4),'L'),2),'p') ; 
       
    end
   end
end

% For all patients
for ii=1:size(ganalysis,1)
    
    % Only process on the first 
    gana = ganalysis{ii,1}; 
    
    % Finds idx in the localization table 
%     idx_subjloc = find(ismember(m_table_as(:,id_subj),upper(gana.subj))) ; 
    % ASD : modif for CLC (2015/12/11)
    idx_subjloc = find(ismember(m_table_as(:,id_subj),gana.subj)) ; 
       
    
%     table_subj = m_table_as(idx_subjloc,:);
%    
    % Filtre sur les patients localized
    if isempty(idx_subjloc); 
        continue;
    end
     
    idx_ganalysisloc = cat(1,idx_ganalysisloc,ii);

    [id, labels,idx] = filter_localized_channels(m_table_as, idx_subjloc , id_contact, id_ncontact, gana.labels) ; 
   
    [un, ia, ic] = unique(labels,'stable');
  
    % Prepare table 
    table_tmp = m_table_as(idx_subjloc(idx),:) ; 
    
    pt_labels_dupli = strcat(m_table_as(idx_subjloc(1),id_subj),'_',labels) ; 
    
    all_labels = cat(1,all_labels,pt_labels_dupli) ; 

%     figure('Name', strcat(gana.subj, 'thresh ' ,num2str(threshp),'dur =',num2str(threshdur)),'Unit','Normalized','Position',[0,0,0.15,1]) ; 
    sig = [] ;
    sigmask = [] ;
    
    table_tmp = cat(2,table_tmp,pt_labels_dupli) ; 
    
    % For all freq.
    for jj=1:size(ganalysis,2)
       
        gana = ganalysis{ii,jj};
        
        tvals = [gana.tvals];
        df = gana.df;
    
        % Get mask of significance 
        h = abs(tvals)>tinv(1-gana.threshp/2,df) ; 
        
        % Removes segments of significance shorter than a threshold duration
        [hf] = filter_timewin_signif(h,gana.threshdur*Fs);
      
        % Get electrodes that show at least one signififcant time point in
        % the post stim peiod
        tmp = sum(hf(:,gana.Time>0)')~=0 ;
        gana.effect = tmp(id(ic));
        
        % Save tvals thresholded
        sig = cat(3,sig,gana.tvals(id(ic),:)) ; 
        sigmask = cat(3,sigmask,gana.tvals(id(ic),:).*hf(id(ic),:)) ; 
       
        tmp =sum(abs(tvals(:,gana.Time>0).*hf(:,gana.Time>0)),2);      
        gana.sumt = tmp(id(ic))'; 

        %% idx : monopolaire  to bipol
        %% ic bipol to mono 
        ganalysis{ii,jj} = gana; 
    
        % Adds bipolar labels in table
        table_tmp = cat(2,table_tmp,arrayfun(@num2cell,gana.effect')); 

    end
    
    s= cat(1,s,sig);
    smask= cat(1,smask,sigmask);
    m_table_effect = cat(1,m_table_effect,table_tmp);

end

