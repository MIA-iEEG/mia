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
% This function browse recursivvely in MarsPower working directory to fill
% out a table containgin all files (processes) that have been produced
function [mia_table, sFiles] = create_table_workdir(varargin)

MAINDIR = varargin{1};
current_loctable = varargin{2};

% Reads all folders that are in MAINDIR 
d = dir(MAINDIR); 
isub = [d(:).isdir]; % returns logical vector if is folder
pt = {d(isub).name}';
pt(ismember(pt,{'.','..'})) = []; % Removes . and ..

ID_PATIENT = 1 ; 
ID_METHOD = 2 ; 
ID_MONTAGE = 3 ; 
ID_FREQS = 4 ; 
ID_SRATE = 5 ; 
ID_REMEVK = 6;
ID_NUMSTATS= 7 ;
ID_CONTACTSLOC = 8 ; 
ID = 9 ;

mia_table = {'','','','','','','',''} ;
sFiles = [] ;
ct = 1;

% Loop through subject from working directrory 
for pp=1:length(pt)

%     % Reads all folders that are in MAINDIR 
	fname = dir(fullfile(MAINDIR,pt{pp},strcat(pt{pp},'*.mat'))); 
  
    % If there is a localisation table for this patient
    if ~isempty(current_loctable) && ~isempty(current_loctable.m_table_all)
%         m_table = load(m_table_fname) ; 
        n_contacts = num2str(sum(strcmp(pt{pp},current_loctable.m_table_all(:,ID_PATIENT)))) ; 
        
    else 
        n_contacts = num2str(0); 
    end
    
    % Loop through the current patient files 
    for ff=1:length(fname) 
        
        % Get all underscores AFTER patient name 
       [~,filename,~]=fileparts(fname(ff).name);
       [mtg remain] = strtok(filename(length(pt{pp})+1:end),'_') ; 
        
       %It is not a montage file
       if isempty(strfind(remain,'_montage'))
        % Valid montage, contnue to analyse
        if strcmpi(mtg,'monopolar') || strcmpi(mtg,'bipolar')
          
            [method remain] = strtok(remain,'_') ;
            
            % Valid methods continue to analyse
            if  strcmpi(method,'LFP') || strcmpi(method,'morlet')
                
                 [datatype remain] = strtok(remain,'_') ;
           
                 % Valid datatypecontinue to analyse
                if  strcmpi(datatype,'data') 
                
                    % Look for a stat file 
                    stat_filename=fullfile(MAINDIR,pt(pp),strrep(fname(ff).name,'_data','_stats'));
                    dat_filename=fullfile(MAINDIR,pt{pp},fname(ff).name) ;
                    
                    % if stat file exist load and count number of stat
                    % computed
                    if exist(char(stat_filename),'file')
                        stat=load(char(stat_filename));
                        mia_table{ct,ID_NUMSTATS}=num2str(length(stat.stats));
                    else
                        mia_table{ct,ID_NUMSTATS}=num2str(0);
                    end
                
                    % Fill out the table (at this point the file is ok to
                    % add in the table)
                    mia_table{ct,ID_PATIENT} =  pt{pp} ;         
                    mia_table{ct,ID_CONTACTSLOC} =  n_contacts ;         
                    mia_table{ct,ID_METHOD} = method;  
                    mia_table{ct,ID_MONTAGE} = mtg;
                    mia_table{ct,ID} = num2str(ct); 

                    if  strcmpi(method,'LFP')  
                        freq = '-';
                        remEvk = '-';
                    else
                        strfreq = strsplit(remain,'_');
                        freq = sprintf('%s-%s Hz (step %s)',strfreq{3},strfreq{4},strfreq{2});
                        remEvk = num2str(double(length(strfreq)==5 && isequal(strfreq{5},'removeEvoked'))); % check if tf decomp processed with removal of evk response or not
                    end
                    mia_table{ct,ID_FREQS} = freq;
                    d = load(dat_filename,'Time') ; 
                    mia_table{ct,ID_SRATE} = sprintf('%d Hz',ceil(1/mean(diff(d.Time))));
                    mia_table{ct,ID_REMEVK} = remEvk;

                    % Fill out the list of file 
                    sFiles{ct} = fullfile(MAINDIR,pt{pp},fname(ff).name);
                    ct = ct +1 ; 
                end
            end
        end
        
        end
    end
        
end 
end
