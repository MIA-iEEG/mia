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
    nUnderscores = numel(bilabels{jj}) - numel(strrep(bilabels{jj},'_',''));

    if nUnderscores >=1

        % If there is only one underscore, we separate the two contacts with the underscore as a delimiter

        if nUnderscores==1
            % Separte bipolar channels labels into mono
            [first,second] = strtok(bilabels(jj),'_') ; 
            second = strrep(second,'_',''); 

        % If there are two underscores, we separate the two contacts with the underscore located in the central position as a delimiter
        elseif nUnderscores==2

            elec_label = bilabels{jj};

            % positions of underscores
            idx_us = strfind(elec_label,'_');

            % first possible split
            left1 = elec_label(1:idx_us(1)-1);
            right1 = elec_label(idx_us(1)+1:end);

            % second possible split
            left2 = elec_label(1:idx_us(2)-1);
            right2 = elec_label(idx_us(2)+1:end);

            % check split 1 against available contacts
            idx1_test_1 = find(ismember(contacts,{left1}));
            idx2_test_1 = find(ismember(contacts,{right1}));

            % check split 2 against available contacts
            idx1_test_2 = find(ismember(contacts,{left2}));
            idx2_test_2 = find(ismember(contacts,{right2}));

            if ~isempty(idx1_test_1) && ~isempty(idx2_test_1)
                first = {left1};
                second = {right1};

            elseif ~isempty(idx1_test_2) && ~isempty(idx2_test_2)
                first = {left2};
                second = {right2};

            elseif inclusive~=0 && (~isempty(idx1_test_1) || ~isempty(idx2_test_1))
                first = {left1};
                second = {right1};

            elseif inclusive~=0 && (~isempty(idx1_test_2) || ~isempty(idx2_test_2))
                first = {left2};
                second = {right2};

            else
            % If neither split matches available contacts, throw an error
                error('error(['MIA error: electrode label contains underscores in an unexpected format. ', ...
                'Please recheck your channel naming. If the problem persists, contact the developers ', ...
                'or open an issue on the MIA GitHub repository.'])')
            end

        % Handle the case where there is more than one underscore 
        elseif nUnderscores==3
            % Get the current bipolar electrode label
            elec_label = bilabels{jj}; 

            % Get the indices of all underscores 
            idx = strfind(elec_label,'_');

            % Separate the strings with the underscore located in the central
            % position
            first = elec_label(1:idx(2)-1) ;
            second = elec_label(idx(2)+1:end) ;

        else
            % Error message to flag the usage of uncommon naming conventions with more than 3 underscores in the bipolar channel label
            error('error(['MIA error: electrode label contains underscores in an unexpected format. ', ...
                'Please recheck your channel naming. If the problem persists, contact the developers ', ...
                'or open an issue on the MIA GitHub repository.'])')
        end
    end

    %Monoplolar case: We just check that the canal has been located in
    %m_table_as
    if isempty(second{1})
        idx1 = find(ismember(contacts,first)) ; 
        if ~isempty(idx1)
            labels = cat(1,labels,bilabels(jj)) ;
            idx = cat(1, idx, idx1);
            id(ct) = jj;
            ct=ct+1;
        end
    else
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