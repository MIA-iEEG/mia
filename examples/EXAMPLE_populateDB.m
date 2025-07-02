function [] = EXAMPLE_populateDB()
% ========================================================================
% This file is part of EXAMPLE project
% 
% Free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This code is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%  
% Copyright (C) 2023 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)
% ------------------------------------------------------------------------
%
% DESCRIPTION : 
%   Function aims to populate a brainstorm database from a directory
%   containging one folder per participant (name of the folder will be used
%   as the name of the subject). (Brainstorm must be up and running at exec)
%   
%   INPUT : Each subject folder must include : 
%       * EEG files : collected with Biosemi system and preprocessed with 
%                    BrainAnalyzer thus three files (.eeg, .vmrk,.vhdr)
%       * <file>_badchan.txt : List of bad channels 
%       * <file>_badtrials.txt : List of bad trials 
%       * <file>_refchan.txt : List of (re)reference channels 
%   OUTPUT : a Brainstorm database 
%
% The processing stages include : 
%   1/ Review raw files (EEG Brainamp : preprocessed files) 
%   2/ Use the Biosemi 64 10-10 electrode cap positions 
%   3/ Reject bad channels (identified previously e.g. in BV) 
%   4/ Import epochs, remove DC offset , resample from 1k to 200Hz
%   5/ Reject bad trials 
%   6/ Compute noise covariance matrix 
%   7/ Compute head model 
%   8/ Compute source reconstruction using dSPM 
%   9/ Compute time-frequency decomposition using Morlet wavelet

%% ===== This block of variables are to be defined by user 
INDIR = '/path/to/dataset'; % Input directory (should contain one folder per participant)
EventsMarkers = {'S 4','S 8'}; % event markers to read in the raw data files
ProtocolName = 'MY_PROTOCOL'; % Name of the new protocol created in Brainstorm

% Parameters
Epoch = [-5, 9.98];     % time window for epoching
New_sfreq = 200;       % new sampling rate 
Baseline_DC = [-5, 9.98];       % baseline substraction ; empty = no DC offset removal
Baseline_tf = [-4, -2]; % time window for TF normalization
Baseline_covar = [-4.999, -1.999];      % time window for noise covariance 
Conductivities = [0.33, 0.004, 0.33];  % conductivities for head model
TimeFreqExplo = 1:80;                  % Frequencies to explore in TF analysis

%% ===== end

% Reads all folders that are in INDIR 
d = dir(INDIR); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

if isempty(subjects)    
    fprintf(sprintf('\n WARNING : NO DATA to process in %s\n',INDIR));
    return;
end


% Starts Brainsorm if it is not already running
if ~brainstorm('status'); brainstorm nogui; end

% Create new protocol
gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);

% Loop through all subjects
for jj=1:length(subjects) 

    % Prepare EEG file name 
    fname= dir(fullfile(INDIR,subjects{jj},strcat(subjects{jj},'.eeg')));
     
    [~, participant_filename,~] = fileparts(fname.name);

    % Check if exist
    if isempty(fname) 
        fprintf(strcat('\nFile not found :\t ',fullfile(subjects{jj},strcat(subjects{jj},'.eeg ------- skip\n'))));
    else

        RawFile = fullfile(INDIR,subjects{jj},fname.name ); 
        
        % Get reference channels file 
        RefChannels = fullfile(INDIR,subjects{jj},strcat(participant_filename,'_refchan.txt')); 
               
        % Get bad channels file 
        BadChannels = fullfile(INDIR,subjects{jj},strcat(participant_filename,'_badchan.txt')); 
        
        % Get bad trials file
        BadTrials = fullfile(INDIR,subjects{jj},strcat(participant_filename,'_badtrials.txt')); 
        
        % Get bad channels per trial
        BadTrialChan = fullfile(INDIR,subjects{jj},strcat(participant_filename,'_badtrialchan.txt')); 
        
        % Call BST functions
        process_pipeline(subjects{jj},RawFile, EventsMarkers, RefChannels, BadChannels, BadTrials, BadTrialChan, Epoch, New_sfreq, Baseline_DC,Baseline_tf, Baseline_covar, Conductivities, TimeFreqExplo)

    end
end
 
%% ===== Function applying the Brainstorm pipeline to the data ===== 
function [] = process_pipeline(SubjectName,RawFile, EventsMarkers, RefChannels, BadChannels, BadTrials, BadTrialChan, Epoch, New_sfreq, Baseline_DC, Baseline_tf, Baseline_covar, Conductivities, TimeFreqExplo)

% Input files
sFiles = [];

% Start a new report
bst_report('Start', sFiles);

% Creates a subject in database using the dfault anatomy and default
% channels
[~, ~] = db_add_subject(SubjectName, [], 1, 0);
panel_protocols('UpdateTree'); % Update the Protocol in GUI

% Process: Create link to raw file
sFiles = bst_process('CallProcess', 'process_import_data_raw', sFiles, [], ...
    'subjectname',    SubjectName, ...
    'datafile',       {RawFile, 'EEG-BRAINAMP'}, ...
    'channelreplace', 1, ...
    'channelalign',   1, ...
    'evtmode',        'value');

% Process: Add EEG positions
sFiles = bst_process('CallProcess', 'process_channel_addloc', sFiles, [], ...
    'channelfile', {'', ''}, ...
    'usedefault',  'ICBM152: BioSemi 64 10-10', ...  % ICBM152: BioSemi 64 10-10
    'fixunits',    1, ...
    'vox2ras',     1, ...
    'mrifile',     {'', ''}, ...
    'fiducials',   []);

% Mark bad channels
if exist(BadChannels,'file') 
     % Process: Set bad channels
    chan_sFiles = bst_process('CallProcess', 'process_channel_setbad', sFiles, [], ...
        'sensortypes', fileread(BadChannels));
else 
    fprintf(sprintf('%s does not exist  ------- No global BAD channels\n',BadChannels));
end

% Mark bad channels
if exist(RefChannels,'file') 
    
    % Get the electrode(s) to re-reference
    ref_chans = fileread(RefChannels); 
    
    if ~isempty(fileread(RefChannels))
        % Re-reference on average sepcified electrodes (if any)
        ref_sFiles = bst_process('CallProcess', 'process_eegref', sFiles, [], ...
            'eegref',      ref_chans, ... % Mastoids
            'sensortypes', 'EEG');
    end
end

% % Save rawFiles in order to import several runs from it
sFiles_raw = sFiles;
sFilesAvg = [];

% Process: Compute head model
sFilesAnat = bst_process('CallProcess', 'process_headmodel', sFiles_raw.FileName, [], ...
    'Comment',     '', ...
    'sourcespace', 1, ...  % Cortex surface
    'meg',         1, ...  
    'eeg',         2, ...  % 3-shell sphere
    'ecog',        1, ...  
    'seeg',        1, ...  
    'openmeeg',    struct(...
         'BemSelect',    [1, 1, 1], ...
         'BemCond',      Conductivities, ...
         'BemNames',     {{'Scalp', 'Skull', 'Brain'}}, ...
         'BemFiles',     {{}}, ...
         'isAdjoint',    0, ...
         'isAdaptative', 1, ...
         'isSplit',      0, ...
         'SplitLength',  4000), ...
         'channelfile', '');
     
% Loop through conditions 
for cc=1:length(EventsMarkers)
    
    % Process: Import MEG/EEG: Events
    sFilesEpochs = bst_process('CallProcess', 'process_import_data_event', ...
        sFiles_raw , [], ...   
        'subjectname',   SubjectName, ...
        'condition',     EventsMarkers{cc}, ...
        'eventname',     EventsMarkers{cc}, ...
        'timewindow',    [], ...
        'epochtime',     Epoch, ...
        'split',         0, ...
        'createcond',    1, ...
        'ignoreshort',   1, ...
        'usectfcomp',    1, ...
        'usessp',        1, ...
        'freq',          New_sfreq, ...
        'baseline',      Baseline_DC, ...
        'blsensortypes', 'MEG, EEG');

    % Get all the runs for this subject (ie the list of the study indices)
    iStudyOther = [sFilesEpochs(1).iStudy];
    % Copy the forward model file to the other runs
    sHeadmodel = bst_get('HeadModelForStudy', sFilesAnat(1).iStudy);
    db_add(iStudyOther, sHeadmodel.FileName);
  
    % Mark bad channels on specific trials
    if exist(BadTrialChan,'file') 

        % Get the electrode(s) to mark as bad in specific trials
        trial_chan = readtable(BadTrialChan); 

        % Get channel file
        channelFile = bst_get('ChannelFileForStudy',  sFiles.FileName); 
        sChannel = in_bst_channel(channelFile);

        fileNames = {sFilesEpochs.FileName}; 
        
        % List of bad channels per file
        for iFile = 1 : length(trial_chan.TrialID)    
    
            % Reads epoch in which there are one or several channels to mark as bad
            sFile = in_bst_data(fileNames{trial_chan.TrialID(iFile)}, 'ChannelFlag', 'Comment');    
           
            % Process: Set bad channels
            chan_sFiles = bst_process('CallProcess', 'process_channel_setbad', sFilesEpochs(trial_chan.TrialID(iFile)), [], ...
            'sensortypes', trial_chan.ChannelList{iFile});
    
        end
    else
        fprintf(sprintf('%s does not exist  ------- No BAD channels marked per trial\n',BadChannels));
    end
    
    % Mark bad trials on all channels
    if exist(BadTrials,'file') 
  
        % Mark bad trials 
        idxBadTrials = str2num(fileread(BadTrials)); 
        SetTrialStatus({sFilesEpochs(idxBadTrials).FileName},1);
    end
    
    % Remove bad trials from subsequent computations
    sFilesEpochs(idxBadTrials) = [];

    % Process: Compute covariance (noise or data)
    sFilesCovar = bst_process('CallProcess', 'process_noisecov', sFilesEpochs, [], ...
        'baseline',       Baseline_covar, ...
        'datatimewindow', [0, 0], ...
        'sensortypes',    'EEG', ...
        'target',         1, ...  % Noise covariance     (covariance over baseline time window)
        'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
        'identity',       0, ...
        'copycond',       0, ...
        'copysubj',       0, ...
        'copymatch',      0, ...
        'replacefile',    1);  % Replace
     
    % Process: Compute sources [2018] and returns all good trials source files
    % (for TF analysis)
    sFilesSources = bst_process('CallProcess', 'process_inverse_2018', sFilesEpochs, [], ...
        'output',  1, ... % Kernel only: shared
        'inverse', struct(...
        'Comment',        'dSPM-unscaled: EEG', ...
        'InverseMethod',  'minnorm', ...
        'InverseMeasure', 'dspm2018', ...
        'SourceOrient',   {{'free'}}, ...
        'Loose',          0.2, ...
        'UseDepth',       1, ...
        'WeightExp',      0.5, ...
        'WeightLimit',    10, ...
        'NoiseMethod',    'reg', ...
        'NoiseReg',       0.1, ...
        'SnrMethod',      'fixed', ...
        'SnrRms',         1e-06, ...
        'SnrFixed',       3, ...
        'ComputeKernel',  1, ...
        'DataTypes',      {{'EEG'}}));

    % Process: Average: By condition : For NOW this is just to keep track
    % of conditions
     avg_sFiles = bst_process(...
        'CallProcess', 'process_average', ...
        sFilesEpochs, [], ...
        'avgtype', 4, ...
        'avg_func', 1, ...  % <HTML>Arithmetic average: <FONT color="#777777">mean(x)</FONT>
        'keepevents', 0);
    
    % Process: Time-frequency (Morlet wavelets)
    sFilesTF = bst_process('CallProcess', 'process_timefreq', sFilesSources, [], ...
        'clusters',      {'Desikan-Killiany', {'postcentral L'}}, ...
        'scoutfunc',     1, ...  % Mean
        'edit',          struct(...
        'Comment',         'Scout postcentral L,Avg,Power,1-80Hz', ...
        'TimeBands',       [], ...
        'Freqs',           TimeFreqExplo, ...
        'MorletFc',        1, ...
        'MorletFwhmTc',    3, ...
        'ClusterFuncTime', 'after', ...
        'Measure',         'power', ...
        'Output',          'average', ...
        'RemoveEvoked',    0, ...
        'SaveKernel',      0), ...
        'normalize2020', 0, ...
        'normalize',     'none');  % None: Save non-standardized time-frequency maps

    % Process: Event-related perturbation (ERS/ERD): [-4.000s,-2.000s]
    sFilesTF = bst_process('CallProcess', 'process_baseline_norm', sFilesTF, [], ...
        'baseline',  Baseline_tf, ...
        'method',    'ersd', ...  % Event-related perturbation (ERS/ERD):    x_std = (x - &mu;) / &mu; * 100
        'overwrite', 0);

     sFilesAvg = [sFilesAvg, avg_sFiles]; 
     
     sFilesEpochsAllCond{cc} = sFilesEpochs; 
     
end

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);
% bst_report('Export', ReportFile, ExportDir);

end
end

%% ===== SET STUDY BAD TRIALS ===== copy from process_detectbad.m (Brainstorm)
% USAGE:  SetTrialStatus(FileNames, isBad)
%         SetTrialStatus(FileName, isBad)
%         SetTrialStatus(BstNodes, isBad)
function SetTrialStatus(FileNames, isBad)
    bst_progress('start', 'Set trial status', 'Updating list of bad trials...');
    % ===== PARSE INPUTS =====
    % CALL: SetTrialStatus(FileName, isBad)
    if ischar(FileNames)
        FileNames = {FileNames};
        [tmp__, iStudies, iDatas] = bst_get('DataFile', FileNames{1});
    % CALL: SetTrialStatus(FileNames, isBad)
    elseif iscell(FileNames)
        % Get studies indices
        iStudies = zeros(size(FileNames));
        iDatas   = zeros(size(FileNames));
        for i = 1:length(FileNames)
            [tmp__, iStudies(i), iDatas(i)] = bst_get('DataFile', FileNames{i});
        end
    % CALL: SetTrialStatus(BstNodes, isBad)
    else
        % Get dependent nodes
        [iStudies, iDatas] = tree_dependencies(FileNames, 'data', [], 1);
        % If an error occurred when looking for the for the files in the database
        if isequal(iStudies, -10)
            bst_error('Error in file selection.', 'Set trial status', 0);
            return;
        end
        % Get study
        sStudies = bst_get('Study', iStudies);
        % Get data filenames
        FileNames = cell(size(iStudies));
        for i = 1:length(iStudies)
            FileNames{i} = sStudies(i).Data(iDatas(i)).FileName;
        end
    end
    
    % Get protocol folders
    ProtocolInfo = bst_get('ProtocolInfo');
    % Get unique list of studies
    uniqueStudies = unique(iStudies);
    % Remove path from all files + Remove all BAD events
    for i = 1:length(FileNames)
        % Remove bad events
        if ~isBad
            DataMat = in_bst_data(FileNames{i}, 'Events');
            isModifiedFile = 0;
            for iEvt = 1:length(DataMat.Events)
                [DataMat.Events(iEvt), isModifiedEvt] = panel_record('SetEventGood', DataMat.Events(iEvt), DataMat.Events);
                if isModifiedEvt
                    isModifiedFile = 1;
                end
            end
            if isModifiedFile
                bst_report('Info', 'process_detectbad', FileNames{i}, 'Event names were modified to remove the tag "bad".');
                disp('BST> Event names were modified to remove the tag "bad".');
                bst_save(file_fullpath(FileNames{i}), DataMat, 'v6', 1);
            end
        end
        % Remove path
        [fPath, fBase, fExt] = bst_fileparts(FileNames{i});
        FileNames{i} = [fBase, fExt];
    end
    
    % ===== CHANGE TRIALS STATUS =====
    % Update each the study
    for i = 1:length(uniqueStudies)
        % === CHANGE STATUS IN DATABASE ===
        % Get files for this study
        iStudy = uniqueStudies(i);
        iFiles = find(iStudy == iStudies);
        % Get study
        sStudy = bst_get('Study', iStudy);
        % Mark trial as bad
        [sStudy.Data(iDatas(iFiles)).BadTrial] = deal(isBad);
        % Update database
        bst_set('Study', iStudy, sStudy);
        
        % === CHANGE NODES STATUS ===
        for iFile = 1:length(iFiles)
            % Get node
            bstNode = panel_protocols('GetNode', [], 'data', iStudy, iDatas(iFiles(iFile)));
            % Update node
            if ~isempty(bstNode)
                bstNode.setModifier(isBad);
            end
        end
        
        % === CHANGE STATUS IN STUDY FILE ===
        % Load study file
        StudyFile = bst_fullfile(ProtocolInfo.STUDIES, sStudy.FileName);
        StudyMat = load(StudyFile);
        % Get previous list of bad trials
        if ~isfield(StudyMat, 'BadTrials') || isempty(StudyMat.BadTrials)
            StudyMat.BadTrials = {};
        end
        % Add bad/good trials to current list
        if isBad
            StudyMat.BadTrials = union(StudyMat.BadTrials, FileNames(iFiles));
        else
            StudyMat.BadTrials = setdiff(StudyMat.BadTrials, FileNames(iFiles));
        end
        % Save list of bad trials in the study file
        bst_save(StudyFile, StudyMat, 'v7');
    end
    % Update tree
    %panel_protocols('UpdateNode', 'Study', uniqueStudies);
    panel_protocols('RepaintTree');
    bst_progress('stop');
end



