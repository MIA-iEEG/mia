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
% ASD 2016/05/19: remove overwrite option (append if does not exist)

function [outname] = mia_s5_compute_stats(varargin)

sInputs = varargin{1} ;
OPTIONS = varargin{2} ;

% % Create progress bar for patients processing
% hwait_pt = waitbar(0,'','Units','Normalized','Name','Statistics');
% 
% % Move waitbar up so it does not overlap with second waitbar
% pos = get(hwait_pt,'Position'); pos(2)= pos(2) +0.1;
% set(hwait_pt,'Position',pos);

fname = sInputs ;
[~,filename,~] = fileparts(fname);

% There is NO data file for LFP
if isempty(strfind(filename,'data'))
    zOPTIONS.overwrite='no';
    % Compute zscore and save file
    fname=mia_s3_zscore(cellstr(sInputs(ii)),zOPTIONS);
    fname=char(fname);
end

% Prepare output name
[PATHSTR,filename,~] = fileparts(fname);
[~,patient,~] = fileparts(char(PATHSTR));

% Create output file name
outname = fullfile(PATHSTR,strcat(strrep(filename,'_data','_stats')));

% % Update progress bar (replace _ by \_ for proper display
% waitbar(ii/length(sInputs),hwait_pt,sprintf('Computing %s...',strrep(patient,'_','\_'))) ;

% Init a counter on sets of stats parameters (all contained in stats
% struct) 
ct =1 ; 

% Check if the file have already created and if overwrite or not
if exist(strcat(outname, '.mat'), 'file')
    % Load existing file
    stat=load(outname);

    % Stats is a field (to keep compatibility with older version)
    if isfield(stat,'stats')
        
        % Check if this set of parameters were not used yet
        if mia_check_stats_exist(stat, OPTIONS) 
            fprintf(sprintf('Stats exist in %s\n', outname)); 
            return 
        else 
            ct = length(stat.stats)+1;
        end
    end
end

fprintf('\nComputing Stats. Patient = %s  %d / %d\n\n',patient , 1, length(sInputs) );

load(fname,'Time','zs','zbaseline','labels','freqb' );

% % Compute t he t-test
% [tvals,pvals] = mia_compute_ttest(zs);
Fs=1/(Time(2)-Time(1));

% Compute the duration threshold for significance (threshdur)
[threshdur, ~, ~] = mia_get_bootthresh(OPTIONS.nboot, zs(:,(Time>OPTIONS.baseline(1))&(Time<=OPTIONS.baseline(2)),:), OPTIONS.alpha) ;

% Saves structure
stat.stats(ct).fname = fname;
stat.stats(ct).pthresh = OPTIONS.alpha;
stat.stats(ct).nboot = OPTIONS.nboot;
stat.stats(ct).baseline= OPTIONS.baseline ;
stat.stats(ct).threshdur= threshdur/Fs ;

save(outname,'-struct','stat');

% % Close progress bar
% delete(hwait_pt) ;

end

function [stats_exist] = mia_check_stats_exist(stat,OPTIONS) 

% Check if that combination pthresh/nboot exist
tmp(1,:) = [stat.stats(:).pthresh] ;
tmp(2,:) = [stat.stats(:).nboot] ;
baselines = [stat.stats(:).baseline] ; 
% All beginings of baseline
tmp(3,:) = baselines(1:2:end) ;
% All ends of baseline
tmp(4,:) = baselines(2:2:end) ;

stats_exist = sum(sum(tmp==repmat([OPTIONS.alpha ;OPTIONS.nboot ; OPTIONS.baseline(1) ; OPTIONS.baseline(2) ],1,size(tmp,2)))==4);
          
end

