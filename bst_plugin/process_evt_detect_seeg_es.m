function varargout = process_evt_detect_seeg_es( varargin )
% PROCESS_EVT_DETECT_SEEG_ES: Epileptic spikes detection for a group of recordings file
%
% USAGE:  OutputFiles = process_evt_detect_seeg_es('Run', sProcess, sInputs)

% @=============================================================================
% This process mainly aims to detect "bad" trials using the criteria from Hirshorn et al. (2016) and Staresina et al. 
% (2012, 2016).
% 
% Copyright (c)2000-2020 University of Southern California & McGill University
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
% Authors: Shuai Wang, Anne-Sophie Dubarry

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Detect epileptic artifacts';
    sProcess.FileTag     = 'ES';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Events';
    sProcess.Index       = 46;  % the same Index as PROCESS_EVT_DETECT_BADSEGMENT
    sProcess.Description = 'TBA.';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Sensor types
    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'SEEG';    
    % apply montage
    sProcess.options.montage.Comment = 'Montage name: ';
    sProcess.options.montage.Type    = 'montage';
    sProcess.options.montage.Value   = '';  
    % montage comment
    sProcess.options.label21.Comment = '<I><FONT color="#777777">Montage files calculated in this process would be deleted at the end.</FONT></I>';
    sProcess.options.label21.Type    = 'label';
    % Separator
    sProcess.options.separator1.Type = 'separator';
    sProcess.options.separator1.Comment = ' ';    
    % Title
    sProcess.options.label11.Comment = '<BR><B>Examine Peak and Gradient:</B>:';
    sProcess.options.label11.Type    = 'label';
    % Threshold0 - compare a peak with the Mean amplitude across trials, in terms of SD (Hirshorn et al., 2016)
    sProcess.options.threshold0.Comment = 'The difference between a Peak and the Mean of peaks across trials is greater than: ';
    sProcess.options.threshold0.Type    = 'value';
    sProcess.options.threshold0.Value   = {5, 'SDs', 0};    
    % Threshold1 - exmaine outliers within trial (Staresina et al., 2012, 2016)
    sProcess.options.label22.Comment = 'At least one time point with value and gradient that are beyond: ';
    sProcess.options.label22.Type    = 'label';    
    sProcess.options.threshold1.Comment = {'None    ', 'Inner Fence (Q1-1.5IQ, Q3+1.5IQ)    ', 'Upper Fence (Q1-3IQ, Q3+3IQ)', ''};
    sProcess.options.threshold1.Type    = 'radio_line';
    sProcess.options.threshold1.Value   = 1;
    % Threshold2 - absolute upper limit
    sProcess.options.threshold2.Comment = 'Absolute Peak is greater than: ';
    sProcess.options.threshold2.Type    = 'value';
    sProcess.options.threshold2.Value   = {350, 'uV', []};    
    % Threshold3 - variation between consecutive time points
    sProcess.options.threshold3.Comment = 'Absolute Gradient is greater than: ';
    sProcess.options.threshold3.Type    = 'value';
    sProcess.options.threshold3.Value   = {25, 'uV', []}; 
    % Separator
    sProcess.options.separator2.Type = 'separator';
    sProcess.options.separator2.Comment = ' ';   
    % Title
    sProcess.options.label12.Comment = '<BR><B>Trial rejection</B>:';
    sProcess.options.label12.Type    = 'label';    
    % Threshold4 - bad channel threshold : If a channel led to more than X% of the trials being rejected, this channel was instead rejected.
    sProcess.options.threshold4.Comment = 'Examine channels that led to the rejetion of trials more than: ';
    sProcess.options.threshold4.Type    = 'value';
    sProcess.options.threshold4.Value   = {50, '%', []};     
    % Validation with visually detected BAD events
    SelectOptions = {...
        '', ...                               % Filename
        '', ...                               % FileFormat
        'open', ...                           % Dialog type: {open,save}
        'Import events...', ...               % Window title
        'ImportData', ...                     % LastUsedDir: {ImportData,ImportChannel,ImportAnat,ExportChannel,ExportData,ExportAnat,ExportProtocol,ExportImage,ExportScript}
        'single', ...                         % Selection mode: {single,multiple}
        'files', ...                          % Selection mode: {files,dirs,files_and_dirs}
        bst_get('FileFilters', 'events'), ... % Get all the available file formats
        'EventsIn'};                          % DefaultFormats: {ChannelIn,DataIn,DipolesIn,EventsIn,MriIn,NoiseCovIn,ResultsIn,SspIn,SurfaceIn,TimefreqIn
    sProcess.options.badevtfile.Comment = 'Select event file with visually marked BAD segments: ';
    sProcess.options.badevtfile.Type    = 'filename';
    sProcess.options.badevtfile.Value   = SelectOptions;  
    % reject trials
    sProcess.options.ismarkbadtrials.Comment = 'Reject trials';
    sProcess.options.ismarkbadtrials.Type    = 'checkbox';
    sProcess.options.ismarkbadtrials.Value   = 0;
    % mark trials
    sProcess.options.isaddevent.Comment = 'Add the reason for rejection as event';
    sProcess.options.isaddevent.Type    = 'checkbox';
    sProcess.options.isaddevent.Value   = 0;    
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    n = length(sInputs);  % number of trials
    % ===== GET OPTIONS ===== 
    % define channels
    SensorTypes = sProcess.options.sensortypes.Value;
    % montage name
    MontageName = sProcess.options.montage.Value;
    % get thresholds
    ThresXSDs = sProcess.options.threshold0.Value{1};         % > mean + X SDs (across-trial)
    FenceType = sProcess.options.threshold1.Value;            % identify outliers within a trial
    ThresPeak = sProcess.options.threshold2.Value{1} / 10e5;  % peak; convert microvolts to volts
    ThresGrad = sProcess.options.threshold3.Value{1} / 10e5;  % gradient; convert microvolts to volts  
    ThresChan = sProcess.options.threshold4.Value{1};         % suspected BAD hannels
    % read BAD_suspected events
    EventFile  = sProcess.options.badevtfile.Value{1};  % get filenames to import
    if exist(EventFile,'file')
      %FileFormat = sProcess.options.badevtfile.Value{2};  % only support .csv for now...
      [EventType, EventTime, EventBadTime] = eventscsv(EventFile);      
    else
      EventBadTime = [];
    end    
    % if mark bad trials or only output reports
    ismarkbadtrials = sProcess.options.ismarkbadtrials.Value;
    isaddevent      = sProcess.options.isaddevent.Value;
    
    % ===== APPLY MONTAGE =====
    fprintf('Calculate temporary montage [%s] for trial rejection.\n', MontageName);
    rawInputs = sInputs;
    MontageFiles    = bst_process('CallProcess', 'process_montage_apply', {sInputs.FileName}, [], 'montage', MontageName, 'createchan', 1);    
    sInputs         = bst_process('GetInputStruct', {MontageFiles.FileName});  % use montage files (default: bipolar2) as inputs    
    fprintf('------------------------------------------------------------------------------------------------------\n');
    
    % Get current progressbar position
    progressPos = bst_progress('get');    
    
    % ===== INITIALIZE VECTORS =====    
    % select good channels - all trials have the same configuration of channels
    ChannelConfig    = in_bst_channel(sInputs(1).ChannelFile);
    ChannelSelected  = channel_find(ChannelConfig.Channel, SensorTypes);
    ChannelNames     = {ChannelConfig.Channel(:).Name}';
    ChannelMarks     = in_bst_data(sInputs(1).FileName, 'ChannelFlag');
    ChannelGood      = ChannelMarks.ChannelFlag == 1;
    ChannelGood      = ChannelSelected(ChannelGood(ChannelSelected));
    ChannelGoodNames = ChannelNames(ChannelGood);
    
    % initialize the trial rejection vectors
    nCh = length(ChannelGood);
    iOk = false(1, n);
    ChannelTrialsPeaks    = zeros(nCh, n);
    ChannelTrialsFenceOut = zeros(nCh, n);
    ChannelTrialsPeaksOut = zeros(nCh, n);
    ChannelTrialsGradsOut = zeros(nCh, n);
    
    % ===== For each trial =====
    for iFile = 1:n
        % Progress bar
        bst_progress('text', 'Reading trial to process...');

        % Load epochs        
        DataES = in_bst_data(sInputs(iFile).FileName, 'F', 'Time', 'ChannelFlag', 'Comment');        
        [DataES.TrialType, DataES.TrialIndex] = commentsplit(DataES.Comment);
        DataES.ChannelGood = ChannelGoodNames;
        DataES.FGood       = DataES.F(ChannelGood, :);        
        DataES.TimeN       = length(DataES.Time);
        % check overlap with the bad segments
        if ~isempty(EventBadTime)
          [DataES.TrialEpoch, DataES.BadOverlaps] = badvalidate(DataES.TrialType, DataES.TrialIndex, DataES.Time, EventType, EventTime, EventBadTime);
        else
          DataES.BadOverlaps = [];
        end
        
        % extracts Peaks
        [DataES.FPeaks, DataES.FPeaksIndex] = max(abs(DataES.FGood), [], 2);
        ChannelTrialsPeaks(:, iFile) = DataES.FPeaks;
        % calculate temporal gradient (variation threshold between consecutive points)
        DataES.FGrad = diff(DataES.FGood, 1, 2);
        % identify outliers by examing fences (cf. https://www.itl.nist.gov/div898/handbook/prc/section1/prc16.htm)
        if FenceType ~= 1
          FenceWeight = 1.5 * FenceType;
          % calculate Q1, Q3 and IQ
          DataES.FGoodQ1 = prctile(DataES.FGood, 25, 2);
          DataES.FGoodQ3 = prctile(DataES.FGood, 75, 2);
          DataES.FGoodIQ = DataES.FGoodQ3 - DataES.FGoodQ1;
          DataES.FGradQ1 = prctile(DataES.FGrad, 25, 2);
          DataES.FGradQ3 = prctile(DataES.FGrad, 75, 2);
          DataES.FGradIQ = DataES.FGradQ3 - DataES.FGradQ1;
          % identify outlier time-points
          Outlier1 = (DataES.FGood < (DataES.FGoodQ1 - DataES.FGoodIQ * FenceWeight)) + (DataES.FGood > (DataES.FGoodQ3 + DataES.FGoodIQ * FenceWeight));
          Outlier2 = (DataES.FGrad < (DataES.FGradQ1 - DataES.FGradIQ * FenceWeight)) + (DataES.FGrad > (DataES.FGradQ3 + DataES.FGradIQ * FenceWeight));
          DataES.checkFenceOuts  = logical(sum(Outlier1(:, 2:end) .* Outlier2, 2));
          DataES.ChannelBadFence = DataES.ChannelGood(DataES.checkFenceOuts);
        else
          DataES.checkFenceOuts  = 0; 
          DataES.ChannelBadFence = 0;
        end
        ChannelTrialsFenceOut(:, iFile) = DataES.checkFenceOuts;
        % abs. peak amplitude greater than the upper limit (350uV by default)
        if ThresPeak > 0
          DataES.checkThresPeak  = DataES.FPeaks > ThresPeak;
          DataES.ChannelBadPeaks = DataES.ChannelGood(DataES.checkThresPeak);
        else
          DataES.checkThresPeak  = 0;
          DataES.ChannelBadPeaks = 0;
        end
        ChannelTrialsPeaksOut(:, iFile) = DataES.checkThresPeak;
        % abs. gradient greater than the upper limit (25uV by default)
        if ThresGrad > 0
          DataES.checkThresGrad  = any(DataES.FGrad > ThresGrad, 2);
          DataES.ChannelBadGrads = DataES.ChannelGood(DataES.checkThresGrad);
        else
          DataES.checkThresGrad  = 0;
          DataES.ChannelBadGrads = 0;
        end
        ChannelTrialsGradsOut(:, iFile) = DataES.checkThresGrad;                          
        
        % Progress bar
        bst_progress('text', 'Saving results...');
        bst_progress('set', progressPos + round(3 * iFile / n / 3 * 100));
        % save DataES to the trial structure
        save(file_fullpath(sInputs(iFile).FileName), 'DataES', '-append');        
    end
    
    % reject trials with peak > x SDs away from the mean across all trials
    ChannelTrialsMean = mean(ChannelTrialsPeaks, 2);
    ChannelTrialsStds = std(ChannelTrialsPeaks, [], 2);
    ChannelTrialsPeaksSDs = ChannelTrialsPeaks > ThresXSDs .* ChannelTrialsStds + ChannelTrialsMean;
    
    fprintf('\n======================================== Trial-wise Report ========================================\n');
    for iTrial = 1:n
      % read ES data
      load(file_fullpath(sInputs(iTrial).FileName), 'DataES');
      % read events of monopolar file
      load(file_fullpath(rawInputs(iTrial).FileName), 'Events');

      DataES.checkPeaksXSDs = ChannelTrialsPeaksSDs(:, iTrial);
      if any(DataES.checkPeaksXSDs)
        DataES.ChannelBadXSDs = DataES.ChannelGood(DataES.checkPeaksXSDs);   
      else
        DataES.ChannelBadXSDs = 0;
      end
      
      % report the detected event and related channels
      if any(DataES.BadOverlaps) || any(DataES.checkFenceOuts) || any(DataES.checkThresPeak) || any(DataES.checkThresGrad) || any(DataES.checkPeaksXSDs)
        % number of channels that have bad trials
        DataES.ChannelBadN = sum((DataES.checkFenceOuts + DataES.checkThresPeak + DataES.checkThresGrad + DataES.checkPeaksXSDs) ~= 0);
        % report detected trial and corresponding channels                       
        criterion = '[Trial rejected]:';
        chan2disp = '';
        if any(DataES.BadOverlaps); criterion = ' it overlaps with the BAD segments,'; end
        if any(DataES.checkFenceOuts)
          criterion = [criterion, ' it has at least one time-point as outlier,']; 
          chan2disp = sprintf('[%s] ', DataES.ChannelBadFence{:});
        end
        if any(DataES.checkThresPeak)
          criterion = [criterion, sprintf(' its peak (abs.) > %d uV,', sProcess.options.threshold2.Value{1})];
          chan2disp = [chan2disp, '|', sprintf('[%s] ', DataES.ChannelBadPeaks{:})]; 
        end
        if any(DataES.checkThresGrad)
          criterion = [criterion, sprintf(' its max gradient (abs.) > %d uV,', sProcess.options.threshold3.Value{1})];
          chan2disp = [chan2disp, '|', sprintf('[%s] ', DataES.ChannelBadGrads{:})];
        end                
        if any(DataES.checkPeaksXSDs)
          criterion = [criterion, ' its peak is an outlier among all trials,']; 
          chan2disp = [chan2disp, '|', sprintf('(%s), ', DataES.ChannelBadXSDs{:})]; 
        end  
        criterion = [criterion, sprintf(' and its artifacts are detected in %d channels', DataES.ChannelBadN)];        
        fprintf('Detect Trial %s (#%d) as bad trial because %s: %s. \n', DataES.TrialType, DataES.TrialIndex, criterion, chan2disp);    
        % add ES as an event marker to this trial
        if isaddevent
          if isempty(Events)
            Events(1).label = criterion;                    % mark the reason to reject this trial
            Events(1).color = [0, 0.5, 0];                  % dark green
            Events(1).epochs = 1;
            Events(1).times = 0;                            % mark at 0 ms
            Events(1).reactTimes = [];
            Events(1).select = 1;
            Events(1).channels = {[]};
            Events(1).notes = {['Channels: ', chan2disp]};  % bad channels
          else
            % remove old comments
            oldEvts = cell2mat(cellfun(@(x) contains(x, 'Trial'), {Events.label}, 'UniformOutput', 0));
            Events(oldEvts) = [];
            % add criterion
            iEvt = length(Events) + 1;
            Events(iEvt).label = criterion;                    % mark the reason to reject this trial
            Events(iEvt).color = [0, 0.5, 0];                  % dark green
            Events(iEvt).epochs = 1;
            Events(iEvt).times = 0;                            % mark at 0 ms
            Events(iEvt).reactTimes = [];
            Events(iEvt).select = 1;
            Events(iEvt).channels = {[]};
            Events(iEvt).notes = {['Channels: ', chan2disp]};  % bad channels            
          end
        end
        save(file_fullpath(rawInputs(iTrial).FileName), 'Events', '-append');  % save events to the monopolar trial
      elseif any(DataES.BadOverlaps)
        fprintf('Trial %s (#%d) is in BAD segments but not detected by our method.\n', DataES.TrialType, DataES.TrialIndex)
      else
        iOk(iTrial) = true;
      end
      % save ES data
      save(file_fullpath(sInputs(iTrial).FileName), 'DataES', '-append');
    end
    
    % report the percent of rejected trials
    ChannelTrialsOutAll =  ChannelTrialsFenceOut + ChannelTrialsPeaksOut + ChannelTrialsGradsOut + ChannelTrialsPeaksSDs;
    TrialsOutN = nnz(sum(ChannelTrialsOutAll));
    fprintf('------------------------------------------------------------------------------------------------------\n');
    fprintf('%d out of %d trials are rejected in total. \n', TrialsOutN, n);
    fprintf('%d trials are rejected because the peak is %d SDs away from the mean across trials. \n', nnz(sum(ChannelTrialsPeaksSDs)), ThresXSDs);
    fprintf('%d trials are rejected because of temporal outliers within trial. \n', nnz(sum(ChannelTrialsFenceOut)));
    fprintf('%d trials are rejected because the peak is higher than %d uV. \n', nnz(sum(ChannelTrialsPeaksOut)), sProcess.options.threshold2.Value{1});
    fprintf('%d trials are rejected because the max gradient is higher than %d uV. \n', nnz(sum(ChannelTrialsGradsOut)), sProcess.options.threshold3.Value{1});
    fprintf('------------------------------------------------------------------------------------------------------\n');
    
    % report bad channels
    fprintf('\n======================================== Channel-wise Report ========================================\n');
    fprintf('Out of %d trials : \n', n);
    chan2reject = [];
    ChannelTrialsOutN = sum(logical(ChannelTrialsOutAll), 2);  % number of rejected trials per channel
    [TrialsSortN, TrialsSortI] = sort(ChannelTrialsOutN, 'descend');
    for i = 1:length(ChannelGoodNames)
      %iChannelTrialsOutAll = ChannelTrialsOutAll(iChannel, :);
      %iChannelTrialsOutN = nnz(iChannelTrialsOutAll);
      iChannel = TrialsSortI(i);
      iChannelTrialsOutN = TrialsSortN(i);
      iChannelTrialsOutRatio = iChannelTrialsOutN / n;
      if iChannelTrialsOutRatio > 0
        fprintf('%.2f%% trials are rejected due to the channel %s. \n', iChannelTrialsOutRatio*100, ChannelGoodNames{iChannel});
        if iChannelTrialsOutRatio >= ThresChan/100
          chan2reject = [chan2reject, iChannel];
        end
      end
    end
    if ~isempty(chan2reject)
      ChannelRejectNames  = ChannelGoodNames(chan2reject);
      ChannelTrialsOutNew = ChannelTrialsOutAll;
      ChannelTrialsOutNew(chan2reject, :) = [];  % remove bad channels
      fprintf('In summary, each of the %d channels (%s) led to at least %d%% trials being rejected. \n', length(chan2reject), sprintf('[%s]', ChannelRejectNames{:}), ThresChan);
      fprintf('After removing the channels listed above, the number of rejected trials should be %d. \n', nnz(sum(ChannelTrialsOutNew)));
    else
      fprintf('No channel led to more than %d%% trials being rejected. \n', ThresChan);
    end
    
    % Return all the input files and bad files
    if ismarkbadtrials
      OutputFiles = {rawInputs(iOk).FileName};
      BadFiles    = {rawInputs(~iOk).FileName};
      [BadPaths, BadFiles] = cellfun(@(x) fileparts(file_fullpath(x)), BadFiles, 'UniformOutput', 0);  % extract filenames that are rejected
      BadFiles = cellfun(@(x) sprintf('%s.mat', x), BadFiles, 'UniformOutput', 0);                     % add file extension - .mat
      OutputPaths = unique(BadPaths);
      for iOutputPath = 1:length(OutputPaths)  % usually for each condition
        iBadTrials = cell2mat(cellfun(@(x) strcmp(x, OutputPaths{iOutputPath}), BadPaths, 'UniformOutput', 0));
        BadTrials = BadFiles(logical(iBadTrials))';
        save(fullfile(OutputPaths{iOutputPath}, 'brainstormstudy.mat'), 'BadTrials', '-append');
      end
    end
    % ===== DELETE TEMPORARY FILES =====
    fprintf('Clean up temporary montage [%s] files. \n', MontageName);
    bst_process('CallProcess', 'process_delete', {sInputs.FileName}, [], 'target', 2);  % 'target' = 1 only delete selected files; 2 delete folder
end

%% sub-functions
function [TrialType, TrialIndex] = commentsplit(Comment)
  % this function is specifically used to extract trial type and index from BST trial Comment
  TrialComments = strsplit(Comment);
  TrialOrder = strsplit(TrialComments{2},{'(','#',')'});
  TrialType = TrialComments{1};         % trial type or condition  
  TrialIndex = str2double(TrialOrder{2});  % within-condition index of this trial
end

function [EventType, EventTime, EventBadTime] = eventscsv(EventFile)
  % this function is specifically used to read BST events and bad segments
  % read from CSV
  fid = fopen(EventFile);
    Events = textscan(fid,'%s%f%f','Delimiter',',');
  fclose(fid);
  EventType = Events{1};  % still cell array
  EventTime = [Events{2}, Events{3}];
  % identify BAD segments
  EventBADs = cell2mat(cellfun(@(x) strcmp(x, 'BAD_suspected'), EventType, 'UniformOutput', 0));
  EventBadTime = EventTime(EventBADs, :);
  EventBadTime(:, 2) = sum(EventBadTime, 2);  % start point and end point
end

function [ThisTrial, Overlaps]=badvalidate(TrialType, TrialIndex, TrialWindow, EventType, EventTime, EventBadTime)
  % to check the overlap between a certain trial and the bad segments
  Trials = cell2mat(cellfun(@(x) strcmp(x, TrialType), EventType, 'UniformOutput', 0));
  TrialsStartTime = EventTime(Trials, 1);
  TrialsTime = [TrialsStartTime + TrialWindow(1), TrialsStartTime + TrialWindow(end)];
  ThisTrial = TrialsTime(TrialIndex, :);
  Overlaps = ~logical((EventBadTime(:,1) > ThisTrial(2)) + (EventBadTime(:,2) < ThisTrial(1)));
end







