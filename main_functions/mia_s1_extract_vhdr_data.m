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
function [OutputFile] = mia_s1_extract_vhdr_data(varargin)

% DATA INPUT INITIALISATION
if nargin<2
    OPTIONS = varargin{1};
    
    % Reads all folders that are in MAINDIR
    d = dir(OPTIONS.maindir);
    isub = [d(:).isdir]; % returns logical vector if is folder
    subjects = {d(isub).name}';
    subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..
    
    for ii=1:length(subjects)
        
        % list all .vhdr files (readable in EEGLAB)
        file = dir(fullfile(OPTIONS.maindir,subjects{ii},strcat('*',OPTIONS.suffix,'.vhdr'))) ;
        
        % ASD To move up out of the loop
        if size(file,1)~=1
            mtg =  sprintf(' %d files found for Patient %s' ,  size(file,1) , subjects{ii} );
            error(strcat('Wrong number of .vhdr files.', mtg));
        end
        
        % List all .vhdr
        sInputs(ii)=cellstr(fullfile(OPTIONS.maindir,subjects{ii},file.name));
    end
    
end

if nargin==2
    
    sInputs=varargin{1};
    OPTIONS=varargin{2};
        
end

%one file import 
if nargin == 3
    sInputs=varargin{1};
    OPTIONS=varargin{2};
    NAME=varargin{3};
end

% Create progress bar during data import
hwait = waitbar(0,'','Units','Normalized');
% Make the waitbar stay on top
set(hwait,'WindowStyle','modal')

subjects_list={};

% Loop through all subjects
for ii=1:length(sInputs)
    
    %get the patient .vhdr file path and name
    [file_PATHSTR,file_NAME,file_EXT] = fileparts(char(sInputs(ii)));
    file_NAME=sprintf('%s%s',file_NAME,file_EXT);
    
    %get the subject name (one file impot)
    if nargin==3
        subject_NAME=char(NAME);
    else
        [~,subject_NAME,~] = fileparts(file_PATHSTR);
    end
    
    % Update progress bar (replace _ by \_ for proper display
    waitbar(ii/length(sInputs),hwait,sprintf('Loading %s... %d / %d',char(strrep(subject_NAME,'_','\_')), ii ,length(sInputs))) ;

    % Prepare output filename
    SDIR = char(fullfile(OPTIONS.outdir,subject_NAME)) ;
    
    % Create Output name
    outname = fullfile(SDIR,strcat(subject_NAME,'_signal_LFP'));
    
    fprintf('\n\nExtracting vhdr data. Patient = %s  %d / %d\n\n',subject_NAME , ii ,length(sInputs) );
    
    % Exist the file or not? Overwrite or not?
    if ((exist (strcat(outname, '.mat'), 'file')) && (strcmpi(OPTIONS.overwrite , 'Yes'))) ||(~exist (strcat(outname, '.mat'), 'file'))
        
        % Load grand avergae data (EEG.data : [nbSamp x nbTrials]
        EEG = pop_loadbv(file_PATHSTR, file_NAME, [], []);
        F = EEG.data ;
        Time = EEG.times/ EEG.srate ;
        
        % Display number of trials and channels
        fprintf('\n\nNumber of trials = %d\n',size(F,3))
        fprintf('Number of channels = %d\n',size(F,1))
        labels = {EEG.chanlocs.labels};
        
        % Removes user specified channels (in bad_chan.xls) 
        if exist(fullfile(OPTIONS.maindir,subject_NAME,'bad_chan.xlsx'),'file')
            % Read bad channels if any
            [~, pdb,~]=xlsread(fullfile(OPTIONS.maindir,subject_NAME,'bad_chan.xlsx'));
            for jj=1:length(pdb)
               % Remove channels
                toRemove = ~cellfun(@isempty,strfind(labels,pdb{jj}));
                labels(toRemove) = [] ;
                F(toRemove,:,:) = [];
            end
            
        end
        
        % Removes user specified channels (in bad_chan.csv file - BIDS compatibility)
        if exist(fullfile(OPTIONS.maindir,subject_NAME,'bad_chan.csv'),'file')
            % Read bad channels if any
            T = readtable(fullfile(OPTIONS.maindir,subject_NAME,'bad_chan.csv'),'ReadVariableNames',0);
            pdb = T.Var1 ; 
            for jj=1:length(pdb)
                % Vire les canaux NULL
                toRemove = ~cellfun(@isempty,strfind(labels,pdb{jj}));
                labels(toRemove) = [] ;
                F(toRemove,:,:) = [];
            end
            
        end
        
        % Systematically removes NULL channels
        isNULL = ~cellfun(@isempty,strfind(labels,'NULL'));
        labels(isNULL) = [] ;
        F(isNULL,:,:) = [];

        % Removes channel for which at least one trial is null  (flat)
        idx_null = find(~min(squeeze(max(abs(F),[],2)),[],2)) ;
        labels(idx_null) = [] ;
        F(idx_null,:,:) = [];
        
        % Reorder labels
        [~,res_index] = mia_sort_nat(labels, 'ascend');
        F=F(res_index,:,:) ; labels = labels(res_index) ;
        Favg = mean(F,3); 
        
           % Create subject folder if does not exist
        if ~exist(SDIR)

           % mkdir(OPTIONS.outdir,subject_NAME) ;
           mkdir(SDIR)

        end
   
        % Save
        save(outname, 'Time', 'F' ,'labels','Favg') ;
        OutputFile{ii} = outname;
        fprintf('\nPatient = %s  %d / %d : Data extracted corectly\n\n',subject_NAME, ii ,length(sInputs) );
        
    end
    
    OutputFile{ii} = outname;
    
end
delete(hwait);

end
