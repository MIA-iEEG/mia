function [id, labels, idx] = mia_filter_localized_channels(m_table_as, idx_subjloc , id_contact, id_ncontact, bilabels, inclusive)
% ***********************************************************************
% ***********************************************************************
%
%  Copyright (C) 2016-2022 CNRS - Universite Aix-Marseille
%
%  This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)
%
% ***********************************************************************

% m_table_as : table of localizations (all patients)
% idx_subjloc : Indices of the current patient in m_table_as
% id_ncontact : Indice of the contacts in m_table_as
% Bilabels : Bipolar channels present in the data
% TODO : ajouter une OPTION pour les canaux bilpolaires : 
% labeling conservateur (seuls les canaux dont les 2 contacts sont dans la
% ROI labélisé ROI
% Labeling restrictif (les canaux bipolaires sont labeliser avec les deux
% ROI de contacts1 et contact2) 

labels = [] ;
ct = 1 ;
idx = [] ;
id  =[] ;

% get contacts from localization table (m_table_as)
% contacts = strrep(strcat(m_table_as(idx_subjloc,id_contact),num2str(cell2mat(m_table_as(idx_subjloc,id_ncontact)))),' ','') ;
contacts = strrep(strcat(m_table_as(idx_subjloc,id_contact),num2str(cell2mat(m_table_as(idx_subjloc,id_ncontact)))),' ','') ;

for jj=1:length(bilabels) 
 
    % Separte bipolar channels labels into mono
    [first,second] = strtok(bilabels(jj),'_') ; 
    second = strrep(second,'_','');
    
    %Cas monoplolaire : On v??rifie juste que le canal a ??t?? localis?? dans
    %m_table_as
    if isempty(second{1})
        idx1 = find(ismember(lower(contacts),lower(first))) ; 
        if ~isempty(idx1)
         labels = cat(1,labels,bilabels(jj)) ;
         idx = cat(1, idx, idx1);
         id(ct) = jj;
         ct=ct+1;
        end
    else
    
        %cas bipolaire : on v??rifie que les deux canaux sont localis??s dans la
        %table et si c'est le cas on d??double la ligne
        % Ne conserve que les indices des canaux bipolaires pour lesquels on a une
        % localisation
        idx1 = find(ismember(lower(contacts),lower(first))) ; 
        idx2 = find(ismember(lower(contacts),lower(second))) ;

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
end
