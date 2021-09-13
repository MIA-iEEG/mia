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
% 2014/3/20 : ASD add input mode in order to be able to run eaither wavelet
% decompo or Hilbert transfos3_compute_TFbandwise(mtg,modetf)
%
% 2016/5/17 : ASD CLEAN code + change output name (no more id)
%
% 2016/5/18  : add progress bar (waitbar) + save history + Times in sec
% everywhere

function [OutputFile] = mia_s4_compute_tf_bandwise(varargin)

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

%Read all the options
if ~isempty(strfind(OPTIONS.modetf,'Morlet')) ; modetf = 'morlet'; end
zbaseline = OPTIONS.zbaseline;
freqb=OPTIONS.freqs  ;
freqstep = freqb(2)-freqb(1);
mtg = lower(OPTIONS.mtg);

% Init output
OutputFile={};

% Create progress bar for patients processing
hwait_pt = waitbar(0,'','Units','Normalized');

% Move waitbar up so it does not overlap with second waitbar
pos = get(hwait_pt,'Position'); pos(2)= pos(2) +0.1;
set(hwait_pt,'Position',pos);

% Loop through all subjects
for ii=1:length(sInputs)
    
    [PATHSTR,NAME,EXT] = fileparts(char(sInputs{ii}));
    
    % Update progress bar (replace _ by \_ for proper display
    islash = strfind(NAME,'_');
    ptname = NAME(1:islash(1)-1) ;
    waitbar(ii/length(sInputs),hwait_pt,sprintf('Computing %s...',strrep(ptname,'_','\_'))) ;
    
    % Prepare output file name
    if OPTIONS.removeEvoked
        output = fullfile(PATHSTR,strrep(NAME,'_signal_LFP', strcat('_',mtg,'_',modetf,'_data_',num2str(freqstep),'_',num2str(freqb(1)),'_',num2str(freqb(end)),'_','removeEvoked'))) ;
    else
        output = fullfile(PATHSTR,strrep(NAME,'_signal_LFP', strcat('_',mtg,'_',modetf,'_data_',num2str(freqstep),'_',num2str(freqb(1)),'_',num2str(freqb(end))))) ;
    end
    
    % Check if the file have already created and if overwrite or not
    if ((exist (strcat(output, '.mat'), 'file')) && (strcmpi(OPTIONS.overwrite , 'Yes'))) ||(~exist (strcat(output, '.mat'), 'file'))
        
        % Load data
        load(sInputs{ii},'Time','F', 'labels');      
        
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
        
        % Sampling rate
        Fs=1/(Time(2)-Time(1));
        
        % Remove evoked response from each trial if selected
        if OPTIONS.removeEvoked
            F = F - mean(F,3);
        end
        
        % Process TF decomposition
        [F,zs] = process_tf(Time,F,freqb,Fs,zbaseline,modetf) ;
        
        % Keep OPTIONS in history
        history = OPTIONS;
        
        % Save results
        save(output,'Time','F','zbaseline','zs','freqb','labels','history');
        
    else
        sprintf('File already exist : %s\n',output) ;
    end
    
    OutputFile= cat(2,OutputFile,strcat(output, '.mat'));
end

delete(hwait_pt) ;
end


function [s, zs] = process_tf(t,dc,freqs,Fs,zbaseline,modetf)

% Time window with no edge
s = zeros(size(dc));
zs = zeros(size(dc));

% Create progress bar
hwait_trials = waitbar(0,sprintf('Computing trial %d/%d...',0, size(dc,3))) ;

% Trial by trial
for trialidx=1:size(dc,3)
    % update progress bar
    waitbar(trialidx/size(dc,3),hwait_trials,sprintf('Computing trial %d/%d...',trialidx, size(dc,3))) ;
    
    if ~strcmpi(modetf,'hilbert')
        [s(:,:,trialidx),zs(:,:,trialidx)] =  compute_wavelet(t,dc(:,:,trialidx),Fs,freqs,zbaseline) ;
    else
        [s(:,:,trialidx),zs(:,:,trialidx)] = compute_hilbert(t,dc(:,:,trialidx),Fs,freqs,zbaseline) ;
    end
end
% Close progress bar
delete(hwait_trials) ;

end


function [s, zs] = compute_wavelet(t,dc,Fs, freqs,zbaseline)

% Loop through contacts 
for contactidx=1:size(dc,1)
    
    % Compute TF decompo
    wt = awt_freqlist(dc(contactidx,:)',Fs,freqs,'Gabor',7);
    
    %Z-score against baseline
    st = abs(wt)';
    baseline = st(:,(t>zbaseline(1))&(t<=zbaseline(2)))' ;
    wtz =  (st - repmat(mean(baseline),length(st),1)')./repmat(std(baseline),length(st),1)';
    
    % Sum abs values of all freq bins for this freq range
    s(contactidx,:)=mean(st);
    zs(contactidx,:)=mean(wtz);
    
    % TODO ASD check that baseline~=0 / for now only display
    if sum(isnan(zs(contactidx,:)))
       fprintf('WARNING : NULL baseline\n');
    end
end
end


function [s, zs] = compute_hilbert(t,dc,Fs,freqs,zbaseline )

% For all contacts
for contactidx=1:size(dc,1)
    % For all freq ranges
    for ff=1:length(freqs)-1
        x = squeeze(dc(contactidx,:)) ;
        [tmp] = bst_bandpass_filtfilt(double(x), Fs, freqs(ff), freqs(ff+1), 0, 'iir') ;
        wt = hilbert(tmp);
        % Normalization : Z-score against baseline
        st(ff,:)= abs(wt);
        baseline = st(ff,(t>zbaseline(1))&(t<=zbaseline(2)))' ;
        wtz (ff,:)=  ( abs(wt) - repmat(mean(baseline),length( abs(wt)),1)')./repmat(std(baseline),length( abs(wt)),1)';
    end
    s(contactidx,:)=mean(st);
    zs(contactidx,:)=mean(wtz);
end
end


