function [id, labels, idx] = mia_filter_localized_channels(m_table_as, idx_subjloc , id_contact, id_ncontact, bilabels, inclusive)
% ***********************************************************************
% ***********************************************************************
%
% Copyright (C) 2016-2022 CNRS - Université Aix-Marseille
%
% This software was developed by
% Anne-Sophie Dubarry (CNRS University of Aix-Marseille)
%
% ***********************************************************************
% m_table_as : table of localizations (all patients)
% idx_subjloc : Indices of the current patient in m_table_as
% id_ncontact : Index of the contacts in m_table_as
% Bilabels : Bipolar channels present in the data
% TODO: add an OPTION for Bilpolar Channels: 
% conservative labeling (only channels where the 2 contacts are in the different ROIs are labeled with both ROI of contacts1 and contact2)
% ROI labelled ROI
% Restrictive labeling (bipolar channels are labeled with both
% ROI of contacts1 and contact2) 
% Contributors : Dewmith Weerasena

labels = [] ;
ct = 1 ;
idx = [] ;
id = [] ;

% get contacts from localization table (m_table_as)
% contacts = strrep(strcat(m_table_as(idx_subjloc,id_contact),num2str(cell2mat(m_table_as(idx_subjloc,id_ncontact)))),' ','') ;
contacts = strrep(strcat(m_table_as(idx_subjloc,id_contact),num2str(cell2mat(m_table_as(idx_subjloc,id_ncontact)))),' ','') ;

for jj=1:length(bilabels) 
   
    % Extract first and second tokens from the bipolar channel label 
    tokens = regexp(bilabels{jj}, '^([^_]+)_([^_]+)', 'tokens');
  
    % Monopolar case but label contains "_"
    if ~isempty(tokens)
        if sum(contains(contacts,tokens{1}{1}))==0
            tokens = []; 
        end
    end

    %Monoplolar case: We just check that the canal has been located in
    %m_table_as
    if isempty(tokens)
        idx1 = find(ismember(contacts,bilabels{jj})) ; 
        if ~isempty(idx1)
            labels = cat(1,labels,bilabels(jj)) ;
            idx = cat(1, idx, idx1);
            id(ct) = jj;
            ct=ct+1;
        end
    else
        first = tokens{1}{1};  % First electrode label (bipolar)
        second = tokens{1}{2}; % Second electrode label (bipolar)

        %bipolar case: we check that the two channels are located in the
        %table and if this is the case we double the line
        % Retains only the indices of bipolar channels for which there is a
        % location
        idx1 = find(ismember(contacts,first)) ; 
        idx2 = find(ismember(contacts,second)) ;

        % At this point if length(idx1) or length(idx2) == 2 - the contact was
        % labelled twice 
        if ~isempty(idx1)&(~isempty(idx2))
            labels = cat(1,labels,bilabels(jj)) ; 
            labels = cat(1,labels,bilabels(jj)) ; 
            idx = cat(1, idx, idx1);
            idx = cat(1, idx, idx2);
            id(ct) = jj;
            ct=ct+1;

        % 2020/01/20 ASD : includes contacts which are present in the data as
        % part of a bilabel channel 
        elseif inclusive~=0 &~isempty(idx1)&isempty(idx2)
            labels = cat(1,labels,bilabels(jj)) ; 
            idx = cat(1, idx, idx1);
            id(ct) = jj;
            ct=ct+1;

        elseif inclusive~=0 &isempty(idx1)&~isempty(idx2)
            labels = cat(1,labels,bilabels(jj)) ; 
            idx = cat(1, idx, idx2);
            id(ct) = jj;
            ct=ct+1;
        end
    end
end