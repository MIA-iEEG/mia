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
% ASD 2016/05/19: remove overwrite option (append if does not exist)

function [OutputFile] = mia_s5_compute_stats(varargin)

if nargin<2
    OPTIONS = varargin{1};
    d = dir(OPTIONS.outdir);
    isub = [d(:).isdir]; % returns logical vector if is folder
    subjects = {d(isub).name}';
    subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..
    for ss=1:length(subjects) ;
        tmp = dir(fullfile(OPTIONS.outdir,subjects{ss},'*_data_*.mat')) ;
        sInputs = cat(2,sInputs,fullfile(OPTIONS.outdir,subjects{ss},{tmp.name})) ;
    end
    
else 
    sInputs = varargin{1} ;
    % If input is only one file, convert to cell for compatibility with
    % multiple data file processing 
    if size(sInputs,1) == 1 ; sInputs = mat2cell(sInputs,[1]) ; end;
    OPTIONS = varargin{2};
end

nboot =OPTIONS.nboot ;
baseline = OPTIONS.baseline ;
OutputFile={};

% Create progress bar for patients processing
hwait_pt = waitbar(0,'','Units','Normalized','Name','Statistics');

% Move waitbar up so it does not overlap with second waitbar
pos = get(hwait_pt,'Position'); pos(2)= pos(2) +0.1;
set(hwait_pt,'Position',pos);

% Loop through all subjects
for ii=1:length(sInputs)
    
    fname = sInputs{ii};
    [~,filename,~] = fileparts(char(fname));
        
    % There is NO data file for LFP
    if isempty(strfind(filename,'data'))
        zOPTIONS.overwrite='no';
        % Compute zscore and save file
        fname=mia_s3_zscore(cellstr(sInputs(ii)),zOPTIONS);
        fname=char(fname);
         
    end
    
    % Prepare output name
    [PATHSTR,filename,EXT1] = fileparts(char(fname));
    [~,patient,EXT] = fileparts(char(PATHSTR));
    
    % Create output file name
    outname = fullfile(PATHSTR,strcat(strrep(filename,'_data','_stats'),EXT));
    
    % Update progress bar (replace _ by \_ for proper display
    waitbar(ii/length(sInputs),hwait_pt,sprintf('Computing %s...',strrep(patient,'_','\_'))) ;
    
    % Check if the file have already created and if overwrite or not
    if ~exist (strcat(outname, '.mat'), 'file')
        ct =1 ;
    else
        % Load existing file
        stat=load(outname);
        
        if isfield(stat,'stats')
            % Check if that combination pthresh/nboot exist
            tmp(1,:) = [stat.stats(:).pthresh] ;
            tmp(2,:) = [stat.stats(:).nboot] ;
            baselines = [stat.stats(:).baseline] ; 
            % All beginings of baseline
            tmp(3,:) = baselines(1:2:end) ;
            % All ends of baseline
            tmp(4,:) = baselines(2:2:end) ;
            if sum(sum(tmp==repmat([OPTIONS.alpha ;OPTIONS.nboot ; OPTIONS.baseline(1) ; OPTIONS.baseline(2) ],1,size(tmp,2)))==4)
                fprintf('SKIP : These parameters already exist\n');
                OutputFile= cat(2,OutputFile,outname);
                continue;
            else
                ct = length(stat.stats)+1;
            end
        else
            ct=1;
        end
    end
    
    fprintf('\nComputing Stats. Patient = %s  %d / %d\n\n',patient , ii, length(sInputs) );
    
    load(fname,'Time','zs','zbaseline','labels','freqb' );
    
    % Compute t he t-test
    [tvals,pvals] = mia_compute_ttest(zs);
    Fs=1/(Time(2)-Time(1));
    
    % Compute the duration threshold for significance (threshdur)
    [threshdur, ~, ~] = mia_get_bootthresh(nboot, zs(:,(Time>baseline(1))&(Time<=baseline(2)),:), OPTIONS.alpha) ;
    
    % Saves structure
    stat.stats(ct).fname = fname;
    stat.stats(ct).pthresh = OPTIONS.alpha;
    stat.stats(ct).nboot = OPTIONS.nboot;
    stat.stats(ct).baseline= OPTIONS.baseline ;
    stat.stats(ct).threshdur= threshdur/Fs ;
    
    save(outname,'-struct','stat');
    
    OutputFile= cat(2,OutputFile,outname);
    
end

% Close progress bar
delete(hwait_pt) ;

end

