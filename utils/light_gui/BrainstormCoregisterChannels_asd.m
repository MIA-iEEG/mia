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
% You should have received a copy of the GNU General Public License
% along with MIA.  If not, see <https://www.gnu.org/licenses/>.
%
% Copyright (C) 2016-2022 CNRS - Université Aix-Marseille
%
% ========================================================================

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
    chanfilename = sStudy(contains({sStudy.Name},'Implantation')).Channel.FileName; % 'Implantation' may not work in all cases and may need adapting

    chandata = load(fullfile(datadir, chanfilename));
    
    % Find current subject's MRI file (expects 1 and only 1):
    mrifileidx = find(cellfun(@(x) contains(x, {'subjectimage_s'})&~contains(x, {'volct'}), {subs.Subject(subidxs(iSubj)).Anatomy.FileName})); 
    % Note that criteria for finding MRI file may change depending on your own naming conventions and the names of files that were imported in BST
    % Alternative looking for MRI file that has been renamed for each subject as SubX_MRI
    % mrifileidx = find(cellfun(@(x) contains(x, {'_MRI'}), {subs.Subject(subidxs(iSubj)).Anatomy.Comment}));
    
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

        % Get MNI coordinates:
        temp(iChan).Loc = cs_convert(mridata, 'scs', 'mni', temp(iChan).Loc')';
        
        % Convert to local coordinates of ICBM152:
        temp(iChan).Loc = cs_convert(coregmridata, 'mni', 'scs', temp(iChan).Loc')';
    
   end
    
    % % Add these channels:
    chanstruct.Channel = [chanstruct.Channel, temp]; %#ok
    
end

% Update COREG channel data:

% Save combined channel file: 
bst_save(fullfile(datadir, 'COREG/@intra/channel.mat'), chanstruct, 'v7');

% Reload all studies:
db_reload_studies(1:bst_get('StudyCount'));

    
