function [m_table_as, status, message] = get_dataloc_table(struct_table, OPTIONS)
% -------------------------------------------------------------------------
% DESCRIPTION
%   Reads localization table
%
% Inputs :
%         filename : Excel file containg the table formatted as follow :

%
% Output:    
%           message : string output message containing doublons (per patient) 
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
% Copyright (C) 2016-2020 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)
% The excel file should contain four columns. Each column should start with
% one of the keywords "Patient", "Electrodes", "Lateralization" and "Region"

m_table_as =[] ;

% Create progress bar for patients processing
hwait_pt = waitbar(0,'','Units','Normalized');

% Move waitbar up so it does not overlap with second waitbar
pos = get(hwait_pt,'Position'); pos(2)= pos(2) +0.1;
set(hwait_pt,'Position',pos);

% Make the waitbar stay on top
set(hwait_pt,'WindowStyle','modal')
status = 1; 
n=1;
for iPt=1:length(struct_table)

    % Update progress bar (replace _ by \_ for proper display
%     waitbar(iPt/length(struct_table),hwait_pt,'Mapping localized contacts with data...');%sprintf('Loading %s...',strrep(u_pt{iPt},'_','\_'))) ;

    pt = struct_table{iPt}.pt ; 
    lat = struct_table{iPt}.lat ; 
    elec = struct_table{iPt}.elec; 
    roi = struct_table{iPt}.roi; 
    
    datafile = fullfile(OPTIONS.maindir,pt,strcat(pt,'_signal_LFP.mat')) ;
    
    if exist(datafile,'file')~=2
        message{n}=sprintf('%s not found : skip.\n',pt);
        status = 0 ;
        n=n+1;
        continue;
        
    else
        % Load labels from data file
%         load(datafile,'labels');
        variableInfo = who('-file', datafile);
        if ismember('isGood',variableInfo)
            load(datafile,'labels','isGood');
            labels = labels(isGood); % Removes contacts that were marked as bad in sanity check
        else
            load(datafile,'labels');
            isGood = ones(1,length(labels));
        end
    end
    %% ASD to fix so that it works for both configuration
    % Preparte a list of electrode labels from the data appending L or R
    % (laterality)
    idx_left_prime = ~logical(cellfun(@isempty,strfind(labels,''''))) ;%subject01 case      
    idx_left_pletter = ~logical(cellfun(@isempty,strfind(labels,'p'))) ; %BRUSE case
    idx_left = idx_left_pletter|idx_left_prime; % for both case of OTp or OT'

    list_data(idx_left) = strcat(labels(idx_left),'L') ;
    list_data(~idx_left) = strcat(labels(~idx_left),'R') ;
    list_data = strrep(list_data,'''','');
    
    % Preparte a list of electrode labels from the data appending L or R
    % (laterality)
    idx_left = logical(strcmp(lat,'L')) ;
    list_table(idx_left) = strcat(elec(idx_left),'L') ;
    list_table(~idx_left) = strcat(elec(~idx_left),'R') ;
    list_table = strrep(list_table,'''','');
    
    % Remove channels that are in the table but not in the data
    bads = ~ismember(lower(list_table),lower(list_data));
%     bads = ~ismember(lower(list_table),lower(list_data(isGood)));
    
    message{n}=sprintf(sprintf('%s : %d contacts found in the data (out of %d listed in the localization table)\n',pt,sum(~bads),length(list_table)));
    n=n+1;
   
    % Display labels of contacts which were found in the table but NOT in
    % the data
    if sum(bads)~=0
        idBads = find(bads) ;
        for bb=1:length(idBads)
            fprintf(sprintf('\n%s_%s found in table but not in data\n',pt,list_table{idBads(bb)}))
        end
    end  
    % if all are bads continue
    if sum(bads)==length(list_table)
        continue;
    end
    
    % Remove bad (channels that are not in the data file)
    elec(bads) = [];
    roi(bads) = [] ;
    lat(bads) = [] ;
    
    % reorder labels (in list of electrodes)
    [~,idx_elec] = sort_nat(lower(elec), 'ascend');
    elec = elec(idx_elec) ;
    roi = roi(idx_elec) ;
    lat = lat(idx_elec) ;
    
    for iCont=1:length(elec);
        matches = regexp(elec(iCont),{'[a-zA-Z]','\d'},'match');
        tab{iCont,1} = pt ;
        tab{iCont,2}= [matches{1}{:}] ;
        tab{iCont,3} = str2num([matches{2}{:}]);
        tab{iCont,4} = cell2mat(lat(iCont)) ;
        tab{iCont,5} = cell2mat(roi(iCont)) ;
        
    end
    
    m_table=tab;
    
    % Load existing file if exist 
    m_table_fname = fullfile(OPTIONS.maindir,pt,'m_table.mat') ; 
    
    % If there is a localisation table for this patient
    if exist(m_table_fname) 
        fprintf(strcat(m_table_fname,' replaced..\n')) ;

    end
    
    % Save patient's table
    save(fullfile(OPTIONS.maindir,pt,'m_table'),'m_table');
    
    % Keep in main table
    m_table_as = cat(1, m_table_as,m_table) ;
    
    clear list_table list_data tab bads
end

% Add 'p' at thend of the contacts labelled if there are left lateralized
endsByP = strcmp(cellfun(@(x)x(end),m_table_as(:,2),'UniformOutput', false),'p') ;
isLeft = strcmp(m_table_as(:,4),'L') ;

idx_needsP = isLeft&~endsByP ; 
m_table_as(idx_needsP,2) = strcat(m_table_as(idx_needsP,2),'p');

if sum(~isLeft&endsByP)~=0 
    fprintf('\nWarning : There are (rigth electrodes) labels ending with ''p'' MIA does not support labels ending with p (confusion possible with left mlateralized channels)');  
end

% % Fix : Replace ' by p
% for iPt=1:size(m_table_as,1)
%     lab = m_table_as(iPt,3) ; 
%     % If it is left ateralized and there is no 'p' already at the end of the label
%     if strcmp(m_table_as(iPt,4),'L') && ~strcmp(lab(end),'p')
%         m_table_as(iPt,2) = strcat(m_table_as(iPt,2),'p');
%     end
% end

%     % Save main table
%     [PATH,~,~]=fileparts(OPTIONS.maindir);
%     save(fullfile(PATH,'m_table_as'),'m_table_as');
% 
    % Remove wait bar
    delete(hwait_pt) ;

end
