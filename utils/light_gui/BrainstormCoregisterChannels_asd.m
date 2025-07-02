
% 
% % List subject in the current BST protocol 
% subjects_bst = bst_get('ProtocolSubjects'); 
% SubjectNames = {subjects_bst.Subject.Name};
% 
% ProtocolInfo = bst_get('ProtocolInfo') ;
% 
% % Process: Select PSD files on each participant
% for iSubj=1:length(SubjectNames)
% 
%          % Process: Select data files in: Subject01/*/Avg: deviant
%             sFiles = bst_process('CallProcess', 'process_select_files_data', [], [], ...
%                     'subjectname',   SubjectNames{iSubj}, ...
%                     'condition',     '', ...
%                     'tag',           '', ...
%                     'includebad',    0, ...
%                     'includeintra',  0, ...
%                     'includecommon', 0);
% 
%         % Load channel file
%         ChannelMat = in_bst_channel(sFiles.ChannelFile);
% 
% end
% 
% 
% 


% CREATE a COREG subject with templatye MNI > ICBM152

%% Channel structure template:
chanstruct = struct('Comment', [], ...  
                    'MegRefCoef', [], ...
                    'Projector', [], ...
                    'TransfMeg', [], ...
                    'TransFMegLabels', [], ...
                    'TransfEeg', [], ...
                    'TransfEegLabels', [], ...
                    'HeadPoints', [], ...
                    'Channel', [], ...
                    'IntraElectrodes', [], ...
                    'History', []);
                
chanstruct.Comment = 'SEEG';

chanstruct.Projector = struct('Comment', [], ...
                              'Components', [], ...
                              'CompMask', [], ...
                              'Status', [], ...
                              'SingVal', []);

chanstruct.HeadPoints = struct('Loc', [], ...
                               'Label', [], ...
                               'Type', []);
chanstruct.History = {};

chanstruct.IntraElectrodes = struct('Name', [], ...
                                    'Type', 'SEEG', ...
                                    'Model', 'DIXI ', ...
                                    'Loc', [], ...
                                    'Color', [0 .8 0], ...
                                    'ContactNumber', 6, ...
                                    'ContactSpacing', .01, ...
                                    'ContactDiameter', .004, ...
                                    'ContactLength', 8e-4, ...
                                    'ElecDiameter', 5e-4, ...
                                    'ElecLength', 0, ...
                                    'Visible', 1);

%% Load and convert all subjects' channels to ICBM152 space:

% Currently loaded protocol:
prot = bst_get('ProtocolInfo'); 

% Data and anatomy directories:
datadir = prot.STUDIES;
anatdir = prot.SUBJECTS;

% All subjects included in the current protocol:
subs = bst_get('ProtocolSubjects');
subnames = {subs.Subject.Name};

% Remove any other COREG...
subs.Subject(~strcmp(subnames, 'COREG') & contains(subnames, 'COREG')) = [];
subnames(~strcmp(subnames, 'COREG') & contains(subnames, 'COREG')) = [];

% Find a subject called 'COREG' (expects this subject to already exist):
coregsubidx = find(cellfun(@(x) strcmp(x, 'COREG'), subnames));
if (isempty(coregsubidx))
    error('Add a new subject to the protol names ''COREG'', with a copy of the ICBM152 MRI and cortex');
end

% Find the MRI file for COREG (expects this to already exist):
mrifileidx = find(cellfun(@(x) contains(x, {'MRI', 'T1'}), {subs.Subject(coregsubidx).Anatomy.FileName}));
if (length(mrifileidx) ~= 1)
    error('Either no or multiple MRI files found for subject %s - cannot proceed', subs.Subject(coregsubidx).Name);
end
mrifilename = subs.Subject(coregsubidx).Anatomy(mrifileidx).FileName;
coregmridata = load(fullfile(anatdir, mrifilename));

% Expand the IntraElectrodes field to the number of existing subjects:
chanstruct.IntraElectrodes = repmat(chanstruct.IntraElectrodes, 1, length(subnames)-1);

% Iterate through each subject, converting channel coordinates. Expects 
% each subject to have 1 and only 1 channel group:
subidxs = setdiff(1:length(subnames), coregsubidx);
chanstruct.Channel = [];

for iSubj = 1:length(subidxs)
    
    % Load current subject's channel structure:
    % Get all studies for this subject
    sStudy = bst_get('StudyWithSubject', subs.Subject(subidxs(iSubj)).FileName) ; 
    chanfilename = sStudy(contains({sStudy.Name},'Implantation')).Channel.FileName;

    % chanfilename = bst_get('StudyWithSubject', subs.Subject(subidxs(iSubj)).FileName).Channel.FileName;

    chandata = load(fullfile(datadir, chanfilename));
    
    % Find current subject's MRI file (expects 1 and only 1):
    mrifileidx = find(cellfun(@(x) contains(x, {'subjectimage_s'})&~contains(x, {'volct'}), {subs.Subject(subidxs(iSubj)).Anatomy.FileName}));

    % Skip if no MRI found
    if isempty(mrifileidx); continue ; end 

    if (length(mrifileidx) ~= 1)
        error('Either no or multiple MRI files found for subject %s - cannot proceed', subs.Subject(subidxs(iSubj)).Name);
    end
    mrifilename = subs.Subject(subidxs(iSubj)).Anatomy(mrifileidx).FileName;
    mridata = load(fullfile(anatdir, mrifilename));

    for k = 1:length(chandata.IntraElectrodes)
        chandata.IntraElectrodes(k).Name = strcat(subs.Subject(subidxs(iSubj)).Name, '_', chandata.IntraElectrodes(k).Name);
    end
    
    % Copy over IntraElectrodes field:
    if iSubj==1
        chanstruct.IntraElectrodes = chandata.IntraElectrodes; % initializae at first patient
    else 
        chanstruct.IntraElectrodes = cat(2,chanstruct.IntraElectrodes,chandata.IntraElectrodes); % concatenates the other ones
    end
    
    % Iterate through each channel:
    temp = chandata.Channel;

    for iChan = 1:length(temp)
        temp(iChan).Comment = subs.Subject(subidxs(iSubj)).Name ;  
        temp(iChan).Name = strcat(subs.Subject(subidxs(iSubj)).Name, '_',temp(iChan).Name);
        temp(iChan).Group= strcat(subs.Subject(subidxs(iSubj)).Name, '_',temp(iChan).Group);

    %     % Set comment to subject's name for later identification:
    %     temp(j).Comment = subs.Subject(subidxs(iSubj)).Name;
    % 
    %     % % Convert to left hemisphere:
    %     % if (temp(j).Loc(2) < 0)
    %     %     temp(j).Loc(2) = -temp(j).Loc(2); 
    %     % end
    %     % 
    %     % Get MNI coordinates:
    %     temp(j).Loc = cs_convert(mridata, 'scs', 'mni', temp(j).Loc')';
    % 
    %     % Convert to local coordinates of ICBM152:
    %     temp(j).Loc = cs_convert(coregmridata, 'mni', 'scs', temp(j).Loc')';

    end
    
    % % Add these channels:
    chanstruct.Channel = [chanstruct.Channel, temp]; %#ok
    
end

% Update COREG channel data:

% Save combined channel file: 
bst_save(fullfile(datadir, 'COREG/@intra/channel.mat'), chanstruct, 'v7');

% Reload all studies:
db_reload_studies(1:bst_get('StudyCount'));

    