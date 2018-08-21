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
% Copyright (C) 2016-2018 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)
% ASD 2016/09/06 : allow the execution of this function with montage in
% argument (directly re-reference if neeeded) + allow input to be OPTIONS
%
% ASD 2017/01/03 : line 84-87 : debug take into account the montage
%
function[OutputFile] = mia_s3_zscore(varargin)

% Check if sInputs is empty and read the files for each patient
if nargin<2
    OPTIONS = varargin{1};
    
    % If subjects are specificed
    if isfield(OPTIONS,'subjects')
        subjects = OPTIONS.subjects;
    else % ALL subjects in directory
        d = dir(OPTIONS.outdir);
        isub = [d(:).isdir]; % returns logical vector if is folder
        subjects = {d(isub).name}';
        subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..
    end
    
    % Prepare input filename
    for ss=1:length(subjects) ;
        tmp = dir(fullfile(OPTIONS.outdir,subjects{ss},'*_signal_LFP.mat')) ;
        sInputs{ss} = fullfile(OPTIONS.outdir,subjects{ss},tmp.name) ;
    end
else % old but we keep for compatibility
    sInputs = varargin{1} ;
    OPTIONS = varargin{2};
end

% Initialize zbaseline (to be saved in each file processed)
zbaseline = OPTIONS.zbaseline ; 

% Create progress bar
hwait = waitbar(0,'','Units','Normalized');

% Loop through all subjects
for ii=1:length(sInputs)
        
    % Prepare output file
    [PATHSTR,NAME,EXT] = fileparts(char(sInputs{ii}));
    output =   fullfile(PATHSTR,strrep(NAME,'_signal_LFP', strcat('_',OPTIONS.mtg,'_LFP_data')));
    
    % Update progress bar (replace _ by \_ for proper display
    waitbar(ii/length(sInputs),hwait,sprintf('Computing %s... %d / %d',char(strrep(NAME,'_','\_')), ii ,length(sInputs))) ;

    % Check if the file have already created and if overwrite or not
    if ((exist (strcat(output, '.mat'), 'file')) && (strcmpi(OPTIONS.overwrite , 'Yes'))) ||(~exist (strcat(output, '.mat'), 'file'))
        
        % Load data
        load(sInputs{ii});

        % Remove Bad channels if needed
        LFP = load(sInputs{ii}) ;         
        if isfield(LFP,'isGood')
            
            F = F(find(LFP.isGood),:,:) ;
            labels=labels(find(LFP.isGood)) ;
        end
        
          % Compute bipolar montage if needed
        if strcmpi(OPTIONS.mtg,'BIPOLAR')
            [F, labels] = mia_make_bipolarmtg(F,labels);
        end
        
        % Calcule zScore
        [t, zs] = process_zs(Time,F,zbaseline);

        % Freq bound = [0 ; Fs]
        freqb = [0 1/(Time(2)-Time(1))] ;
        
        % Save results
        save(output,'Time','F','zbaseline','zs','freqb','labels' );
        
    end
    OutputFile{ii} = strcat(output, '.mat');
end

delete(hwait);

end

function [t, zs] = process_zs(t,s,zbaseline)

% Create progress bar
hwait_trials = waitbar(0,sprintf('Computing trial %d/%d...',0, size(s,3))) ;

% Trial by trial
for trialidx=1:size(s,3)
    
    % update progress bar
    waitbar(trialidx/size(s,3),hwait_trials,sprintf('Computing trial %d/%d...',trialidx, size(s,3))) ;

    % zscore correction
    st = s(:,:,trialidx);
    baseline = st(:,(t>zbaseline(1))&(t<=zbaseline(2)))' ;
    zs(:,:,trialidx) =  (st - repmat(mean(baseline),length(st),1)')./repmat(std(baseline),length(st),1)';
    
end
% Close progress bar
delete(hwait_trials) ;

end


