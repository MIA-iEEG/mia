function [id, labels, idx] = filter_localized_channels(m_table_as, idx_subjloc , id_contact, id_ncontact, bilabels)
% ***********************************************************************
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

% m_table_as : table of localizations (all patients)
% idx_subjloc : Indices of the current patient in m_table_as
% id_ncontact : Indice of the contacts in m_table_as
% Bilabels : Bipolar channels present in the data

labels = [] ;
ct = 1 ;
idx = [] ;
id  =[] ;

% get contacts from localization table (m_table_as)
% contacts = strrep(strcat(m_table_as(idx_subjloc,id_contact),num2str(cell2mat(m_table_as(idx_subjloc,id_ncontact)))),' ','') ;
contacts = strrep(strcat(m_table_as(idx_subjloc,id_contact),num2str(cell2mat(m_table_as(idx_subjloc,id_ncontact)))),' ','') ;

for jj=1:length(bilabels) 

    % Separe les canaux bipolaires en labels mono
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
        
    end
    end

end
