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
% 2014: Modified by JCM to work inside of brainstorm 30-Sep-2015
%
% 2017/05/22 : ASD: Debug (SUBARRY was not the same dimension as DATA which
% caused a crash : replace the indexing way with a boolean 
function [OutputFile] = mia_s1_extract_bst_data(varargin)

if nargin==1
    
    OPTIONS=varargin{1};
    % Reads all folders that are in MAINDIR
    d = dir(OPTIONS.maindir);
    isub = [d(:).isdir]; % returns logical vector if is folder
    subjects = {d(isub).name}';
    subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

    if isempty(subjects)
        error(sprintf('Not data files found in %s',OPTIONS.maindir));
    end
    
else
    OPTIONS=varargin{1};
    subjects=varargin{2};
end

% Create progress bar during data import
hwait = waitbar(0,'','Units','Normalized');
% Make the waitbar stay on top
set(hwait,'WindowStyle','modal')

% Loop through all subjects
for ii=1:length(subjects)
    
    % Prepare output filename
    SDIR = char(fullfile(OPTIONS.outdir,subjects{ii})) ;
    
    % Create subject folder if does not exist
    if ~exist(SDIR)
        mkdir(OPTIONS.outdir,subjects{ii}) ;
    end
    
    % Create Output name
    outname = fullfile(SDIR,strcat(subjects{ii},'_signal_LFP'));
    % Exist the file or not? Overwrite or not?
    if ((exist (strcat(outname, '.mat'), 'file')) && (strcmpi(OPTIONS.overwrite , 'Yes'))) ||(~exist (strcat(outname, '.mat'), 'file'))
                   
        % Update progress bar (replace _ by \_ for proper display
        waitbar(ii/length(subjects),hwait,sprintf('Loading %s... %d / %d',char(strrep(subjects{ii},'_','\_')), ii ,length(subjects))) ;

        % If there is one file to processe
        if nargin==1
            THISDIR = fullfile(OPTIONS.maindir,subjects{ii});
        else % If there is MORE than one file to processe
            THISDIR = fullfile(OPTIONS.maindir);
        end 
        
        CH_FILES = dir(fullfile(THISDIR,'*channel*.mat'));
        DATA_FILES = dir(fullfile(THISDIR,'data*trial*.mat'));
        
        % Load channel file
        CH = load(fullfile(THISDIR,CH_FILES(1).name));
    
        % Load data
        DATA = load(fullfile(THISDIR,DATA_FILES(1).name));
        Time = DATA.Time; % CBB: assume same for all trials

        % Find the valid Subarray channels : Modif ASD 2017/05/23
        bool_subarray = strcmp(OPTIONS.SensorType,{CH.Channel.Type})& (DATA.ChannelFlag==1)';
    
        % Get their corresponding labels
        labels = {CH.Channel(bool_subarray).Name};
         
        % Get the corresponding three dimensional data
        % ASD only keep good channels
        tmp = DATA.F(bool_subarray,:);
        
        F = zeros(size(tmp,1),size(tmp,2),length(DATA_FILES));
        F(:,:,1) = tmp; % first data set
        
        % now the rest
        for jj = 2:length(DATA_FILES)
            DATA = load(fullfile(THISDIR,DATA_FILES(jj).name));
            F(:,:,jj) = DATA.F(bool_subarray,:);
        end
        
        % JCM So now we have the data, the time, and the labels.
        
        % Removes user specified channels (in bad_chan.xls file)
        if exist(fullfile(OPTIONS.maindir,subjects{ii},'bad_chan.xlsx'),'file')
            % Read bad channels if any
            [~, pdb,~]=xlsread(fullfile(OPTIONS.maindir,subjects{ii},'bad_chan.xlsx'));
            for jj=1:length(pdb)
                % Vire les canaux NULL
                isNULL = ~cellfun(@isempty,strfind(labels,pdb{jj}));
                labels(isNULL) = [] ;
                F(isNULL,:,:) = [];
            end
        else
            % Removes "NULL" channels
            isNULL = ~cellfun(@isempty,strfind(labels,'NULL'));
            labels(isNULL) = [] ;
            F(isNULL,:,:) = [];
        end
        
        % Remove channel for which at least one trial is null 
        idx_null = find(~min(squeeze(max(abs(F),[],2)),[],2)) ;
        labels(idx_null) = [] ;
        F(idx_null,:,:) = [];
        
        % Reorder labels
        [~,res_index] = sort_nat(labels, 'ascend');
        F=F(res_index,:,:) ; labels = labels(res_index) ;
        Favg = mean(F,3); 
        
        % Save
        save(outname, 'Time', 'F' ,'labels','Favg') ;
        OutputFile{ii} = outname;
        fprintf('\nPatient = %s  %d / %d : Data extracted corectly\n\n',subjects{ii} , ii ,length(subjects) );
        
    end
    OutputFile{ii} = outname;
end

delete(hwait);

