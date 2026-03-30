function varargout = process_concatenate_channels( varargin )
% PROCESS_SIMULATE_MATRIX: Simulate source signals and saves them as a matrix file.
%
% USAGE:   OutputFiles = process_simulate_sources('Run', sProcess, sInputA)
%               signal = process_simulate_sources('Compute', fnesting, fnested, duration, sRate, couplingPhase, DutyCycle)
 
% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Guiomar Niso, Francois Tadel, 2013-2014

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'MIA: Concatenate Channels';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Standardize';
    sProcess.Index       = 306;
    sProcess.Description = '';

    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'import'};
    sProcess.OutputTypes = {'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 0;

    % === SUBJECT NAME
    sProcess.options.subjectname.Comment = 'New subject name:';
    sProcess.options.subjectname.Type    = 'text';
    sProcess.options.subjectname.Value   = 'COREG';

    % === SUBJECT TO SKIP
    sProcess.options.subskip.Comment = 'Subjects to skip:';
    sProcess.options.subskip.Type    = 'text';
    sProcess.options.subskip.Value   = '';

    % === HELP TEXT AT THE BOTTOM
    sProcess.options.label1.Comment = [ ...
        '<BR>Copies the channels from multiple subject files into a single matrix and creates a new subject.<BR>' ...
        'Subjects which need to be skipped can be listed separated by commas.' ...
    ];
    sProcess.options.label1.Type = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputA) %#ok<DEFNU>
    OutputFiles = {};
    % ===== GET OPTIONS =====
    % Get subject name
    NewSubjectName = file_standardize(sProcess.options.subjectname.Value);
    if isempty(NewSubjectName) 
        bst_report('Error', sProcess, sInputs, 'Subject name is empty.');
        return
    end

    
    % Get protocol name
    subs = bst_get('ProtocolSubjects');
    subskip= strtrim(sProcess.options.subskip.Value);

    % Get the names of all the subjects into an array
    subnames = {subs.Subject.Name};

    % Handle the case of the new subject already exist --> raise an
    % error

    if ismember(NewSubjectName,subnames)
        error('Subject "%s" already exists. Please choose a different name.', NewSubjectName);
    end

    % Add subskip variable with subject names which startes with 'COREG'
    subskip = strjoin([(subnames(startsWith(subnames, 'COREG'))), subskip], ', ');

    % Make sure teh skipped subject names are not repeated
    subskip = unique(subskip, 'stable');


    % Initialize channel structure
    chanstruct.Channel = [];
    
    % Currently loaded protocol:
    prot = bst_get('ProtocolInfo'); 
    
    % Data and anatomy directories:
    datadir = prot.STUDIES;
    anatdir = prot.SUBJECTS;

    
    % coregsubidx = find(cellfun(@(x) strcmp(x, 'COREG'), subnames));
    % if isempty(coregsubidx)
    %     error('Add a new subject to the protol names ''COREG'', with a copy of the ICBM152 MRI and cortex');
    % end
    % 
    % % Find the MRI file for COREG (expects this to already exist):
    % mrifileidx = find(cellfun(@(x) contains(x, {'MRI', 'T1'}), {subs.Subject(coregsubidx).Anatomy.FileName}));
    % if (length(mrifileidx) ~= 1)
    %     error('Please add a subject named ''%s'' with ICBM152 MRI and cortex.', COREG_SUBJECT_NAME);
    % end
    
    % Load COREG MRI once
    %coregmridata = load(fullfile(anatdir, subs.Subject(coregsubidx).Anatomy(mrifileidx).FileName));
    
    subidxs = find(~ismember(subnames,strtrim(strsplit(subskip, ',')))) ; 
    
    % Here check that the subject the user wants to create does not existat
    % already 




    % Then for all other subjects (excpetect the ones we skip)
    for iSubj = 1:length(subidxs)
        
        % Load current subject's Implantation channel 
        sStudy = bst_get('StudyWithSubject', subs.Subject(subidxs(iSubj)).FileName) ; 
        if isempty(sStudy) ; fprintf(strcat('Skipping : ', subs.Subject(subidxs(iSubj)).FileName,'\n')) ; continue ; end 
        chanfilename = sStudy(1).Channel.FileName; 
        %chanfilename = sStudy(find(contains({sStudy.Name},'rawsub'),1)).Channel.FileName; 
        chandata = load(fullfile(datadir, chanfilename));
        
        % % Find current subject's MRI file (expects 1 and only 1):
        % mrifileidx = find(cellfun(@(x) contains(x, {'subjectimage_s'})&~contains(x, {'volct'}), {subs.Subject(subidxs(iSubj)).Anatomy.FileName})); 
        % % Note that criteria for finding MRI file may change depending on your own naming conventions and the names of files that were imported in BST
        % % Alternative looking for MRI file that has been renamed for each subject as SubX_MRI
        % % mrifileidx = find(cellfun(@(x) contains(x, {'_MRI'}), {subs.Subject(subidxs(iSubj)).Anatomy.Comment}));
        % 
        % % Skip if no MRI found
        % if isempty(mrifileidx); continue ; end 
        % 
        % if (length(mrifileidx) ~= 1)
        %     error('Either no or multiple MRI files found for subject %s - cannot proceed', subs.Subject(subidxs(iSubj)).Name);
        % end
        % 
        % mrifilename = subs.Subject(subidxs(iSubj)).Anatomy(mrifileidx).FileName;
        % mridata = load(fullfile(anatdir, mrifilename));
        % 
    
        % Prefix electrode names with subject name
        for k = 1:length(chandata.IntraElectrodes)
            chandata.IntraElectrodes(k).Name = strcat(subs.Subject(subidxs(iSubj)).Name, '_', chandata.IntraElectrodes(k).Name);
        end
    
        % Accumulate IntraElectrodes
        if iSubj==1
            chanstruct = chandata ; chanstruct.Channel = [] ; 
            % chanstruct.IntraElectrodes = chandata.IntraElectrodes; % initializae at first patient
        else 
            chanstruct.IntraElectrodes = cat(2,chanstruct.IntraElectrodes,chandata.IntraElectrodes); % concatenates the other ones
        end
        
        % Iterate through each channel:
        temp =chandata.Channel(strcmp({chandata.Channel.Type},'SEEG')) ; 

        for iChan = 1:length(temp)
            temp(iChan).Comment = subs.Subject(subidxs(iSubj)).Name ;  
            temp(iChan).Name = strcat(subs.Subject(subidxs(iSubj)).Name, '_',temp(iChan).Name) ;
            temp(iChan).Group= subs.Subject(subidxs(iSubj)).Name ;
            % 
            % % Get MNI coordinates:
            % temp(iChan).Loc = cs_convert(mridata, 'scs', 'mni', temp(iChan).Loc')';
            % 
            % % Convert to local coordinates of ICBM152:
            % temp(iChan).Loc = cs_convert(coregmridata, 'mni', 'scs', temp(iChan).Loc')';
            % 
       end
        
        % % Add these channels:
        chanstruct.Channel = [chanstruct.Channel, temp]; 
        
    end

chanstruct.Comment = sprintf('Grand Subject (%d)',size(chanstruct.Channel,2));

[sSubject, iSubject] = db_add_subject(NewSubjectName, [], 0, 0) ; 

% 
% % Color code by patient 
% %% Color Code by Patient
% % Extract a vector of patient IDs (digits following the letter "S")
% vector = cellfun(@(x) sscanf(x, 'S%d'), {chan.Channel.Comment}, 'UniformOutput', true);
% % Get unique patient IDs and map them to sequential integers
% [unique_vals, ~, idx] = unique(vector);
% % Replace all data in dat. F with a patient-based color code
% that. F = repmat(idx, 1, size(that. F, 2));
% % Add a comment in dat
% dat.Comment = sprintf('Patients. Number of patients = %d (Custom colormap)', length(unique_vals));

% Process: Simulate generic signals
sFiles = bst_process('CallProcess', 'process_simulate_matrix', [], [], ...
    'subjectname', NewSubjectName, ...
    'condition',   '', ...
    'samples',     10, ...
    'srate',       1000, ...
    'matlab',      [sprintf('Data(1,:) = sin(2*pi*t); \nData(%d,:) = cos(pi*t) + 1;',size(chanstruct.Channel,2))]);

OutputFiles = import_raw(bst_fullfile(datadir,sFiles.FileName),'BST-MATRIX',iSubject);

% Save combined channel file: 
%bst_save(fullfile(datadir, NewSubjectName, '@intra/channel.mat'), chanstruct, 'v7');
bst_save(fullfile(fileparts(OutputFiles{1}),'channel.mat'), chanstruct, 'v7');

panel_protocols('UpdateTree');

end