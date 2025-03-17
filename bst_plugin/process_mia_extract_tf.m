function varargout = process_mia_extract_tf( varargin )
% PROCESS_RESAMPLE: MIA TF decomposition (by bands)
%
% USAGE:        sInput = process_mia_extract_tf('Run', sProcess, sInput)
% @=============================================================================
% This function was adapted from the Brainstorm software:
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
% Authors: Francois Tadel, A.-Sophie Dubarry 2010-2025

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'MIA: Time-frequency (Morlet by band + 1/f norm)';
    sProcess.FileTag     = 'mia';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Frequency';
    sProcess.Index       = 506;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'results', 'matrix'};
    sProcess.OutputTypes = {'data', 'results', 'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Default values for some options
    sProcess.processDim  = 1;    % Process channel by channel
    sProcess.isSeparator = 1;
    
    % === Freq band
    sProcess.options.freqb.Comment = 'Frequency bands:';
    sProcess.options.freqb.Type    = 'text';
    sProcess.options.freqb.Value   = '50:10:170';
    
    % === Number of cycle
    sProcess.options.ncycle.Comment = 'Number of cycles:';
    sProcess.options.ncycle.Type    = 'text';
    sProcess.options.ncycle.Value   = '7';

    % === Process description
    sProcess.options.label1.Comment = ['Warning: Edges should be removed from baseline <BR>' ...
                                       'Wavelet length <BR>'...
                                       '&nbsp; <B>&sigma;</B> = &nbsp;&nbsp;&nbsp;<FONT color=#7F7F7F>[std(x(iBaseline))]</FONT><BR><BR>'];
    sProcess.options.label1.Type = 'label';
    
    % === Baseline time window
    sProcess.options.baseline.Comment = 'Baseline:';
    sProcess.options.baseline.Type    = 'baseline';
    sProcess.options.baseline.Value   = [];
    sProcess.options.baseline.Group   = 'input';
      
end
  

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sprintf('MIA: %sHz', sProcess.options.freqb.Value);
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    
    % Initialize returned list of files
    OutputFiles = {};

    % Initialize progress bar
    startValue = bst_progress('get');

    % ===== GET INPUT FILES =====
    FileNames = file_fullpath({sInputs.FileName}); %si va

    % Loop over all the input files
    for iInput = 1:length(FileNames) 
        bst_progress('set', round(startValue + (iInput-1) / length(FileNames) * 100));
       
        % ===== OUTPUT STRUCTURE =====
        % Load full original file
        sMat = load(FileNames{iInput});
        sBands = eval(sProcess.options.freqb.Value);
        sFreq = 1 ./ (sMat.Time(2) - sMat.Time(1));

        [TF,zs] = compute_miawavelet(sMat.Time,sMat.F,sFreq,sBands,sProcess.options.baseline.Value{1}, str2num(sProcess.options.ncycle.Value)) ;   
        sMat.F = zs;
        sMat.Comment = [sMat.Comment ' ' '| MIA'];
        fileTag = '_mia';
        
        % Add history entry
        sMat = bst_history('add', sMat, 'process', [func2str(sProcess.Function) ': -------' FormatComment(sProcess) ' | Baseline :' num2str(sProcess.options.baseline.Value{1}) ' | nCycle : ' sProcess.options.ncycle.Value ]);
        
        % ===== SAVE FILE =====
        % Get study description
        sStudy = bst_get('Study', sInputs(iInput).iStudy);
        % Output filename
        OutputFiles{iInput} = file_unique(strrep(FileNames{iInput}, '.mat', [fileTag '.mat']));
        % Save on disk
        bst_save(OutputFiles{iInput}, sMat, 'v6');
        % Register in database
        db_add_data(sInputs(iInput).iStudy, OutputFiles{iInput}, sMat);
    
    end
end


function [s, zs] = compute_miawavelet(vTime, F, sFreq, sBands, zbaseline, nCycle) 
   
    % Initialize results
    s = zeros(size(F));
    zs = zeros(size(F));
  
    for contactidx=1:size(F,1)
        % Compute TF decompo
        wt = mia_awt_freqlist(F(contactidx,:)',sFreq, sBands,'Gabor',nCycle);
        %Z-score against baseline 
        st = abs(wt)'; 
        baseline = st(:,(vTime>zbaseline(1))&(vTime<=zbaseline(2)))' ;
        wtz =  (st - repmat(mean(baseline),length(st),1)')./repmat(std(baseline),length(st),1)';
        % sum abs values of all freq bins for this freq range 
        s(contactidx,:)=sum(st);
        zs(contactidx,:)=sum(wtz);
       
    end
end
  