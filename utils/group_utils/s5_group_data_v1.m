function [ganalysis] = s5_group_data_v1(sFiles,OPTIONS)
% -------------------------------------------------------------------------
% DESCRIPTION
%   Group multiple patients data into one structure
%
% Inputs :
%         grpOPTIONS :
%
% Output:   structure array
%
% Example : ganalysis{1}
%   subj: 'Subject01'
%         fname: 'Subject01_signal_LFP_bipolar_morlet_data_2_10_80_150'
%         tvals: [161x1801 double]
%         pvals: [161x1801 double]
%            df: 240
%         freqb: [80 90 100 110 120 130 140 150]
%          Time: [1x1801 double]
%        labels: {1x161 cell}
%     threshdur: 28
%     dur_array: [1000x1 double]
%       captmax: [1x1000 double]
%             t: [1x1801 double]
%
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
%


ganalysis=[];

% Create a wait bar 
hwait=waitbar(0,sprintf('Processing %d/%d',1,length(sFiles)),'Name','Group analysis');

% Make the waitbar stay on top
set(hwait,'WindowStyle','modal')

% This indice is here because previous versions of MIA was handling several
% col
ii=1;

% for pp=1:length(subj)
for pp=1:length(sFiles)
    
    [~,sSubj,~] = fileparts(fileparts(sFiles{pp})) ; 

    waitbar(pp/length(sFiles),hwait,sprintf('%s %d/%d',strrep(sSubj,'_',' '),pp,length(sFiles)))
    
    % Load data
    data = load(sFiles{pp});
        
    ganalysis{pp,ii}.fname=sFiles{ii};
    ganalysis{pp,ii}.subj=sSubj;
    ganalysis{pp,ii}.df=size(data.F,3);
    ganalysis{pp,ii}.freqb=data.freqb;
    ganalysis{pp,ii}.edges=[OPTIONS.low_twin, OPTIONS.up_twin];
    ganalysis{pp,ii}.threshp=OPTIONS.alpha;
      
    % Create the 'stats' filename
    statfile = strrep(sFiles{ii},'_data','_stats');

    % computation to find or calculate threshdur
    if exist(statfile) %check if stat file already exist

        %load stats file
        stat_data=load(statfile); 

        if isfield(stat_data,'stats') %check for the existance of the field 'stats' (case of 'stat' file not computed)
            for ss=1:length(stat_data.stats)
                if ~isempty(stat_data.stats(ss).fname)&& ~isempty(stat_data.stats(ss).pthresh)&& ~isempty(stat_data.stats(ss).nboot)&& ~isempty(stat_data.stats(ss).threshdur)
                    % If the specific combination of nboot + alpha was already computed 
                    if stat_data.stats(ss).nboot==OPTIONS.nboot && stat_data.stats(ss).pthresh==OPTIONS.alpha 
                        ganalysis{pp,ii}.threshdur=stat_data.stats(ss).threshdur;

                    end
                end
            end
        end
        %No stats has already been calculated with those parameters
        %if no combination match or no stat field
        if ~isfield(ganalysis{pp},'threshdur')
            statfile=mia_s5_compute_stats(cellstr(sFile),OPTIONS); %stats computation
            stat_data=load(char(statfile));
            for ss=1:length(stat_data.stats)
                if ~isempty(stat_data.stats(ss).fname)&& ~isempty(stat_data.stats(ss).pthresh)&& ~isempty(stat_data.stats(ss).nboot)&& ~isempty(stat_data.stats(ss).threshdur)
                    if stat_data.stats(ss).nboot==OPTIONS.nboot && stat_data.stats(ss).pthresh==OPTIONS.alpha
                        ganalysis{pp,ii}.threshdur=stat_data.stats(ss).threshdur; %keep the new calculated threshdur

                    end
                end
            end
        end

    else  %cas of inexistent statfile
        statfile=mia_s5_compute_stats(sFiles{ii},OPTIONS);
        stat_data=load(char(statfile));
        ganalysis{pp,ii}.threshdur=stat_data.stats.threshdur;
    end

    % Compute Ttest
    [tvals,pvals] = mia_compute_ttest(data.zs);

    % Remove value of Time, tvals and pvals that are not in the time
    %window
    for ss=1:length(data.Time)
        if data.Time(ss)<=OPTIONS.low_twin || data.Time(ss)>=OPTIONS.up_twin
            data.Time(ss)=0;
        end
    end
    % Keep pvals, tvals and Time datas without edges effect
    ganalysis{pp,ii}.Time=data.Time(find(data.Time));
    ganalysis{pp,ii}.tvals=tvals(:,find(data.Time));
    ganalysis{pp,ii}.pvals=pvals(:,find(data.Time));

    % Place the p(for left contacts) at the right place on contact name
    data.labels=strrep(data.labels,'''','p');
    ganalysis{pp,ii}.labels=data.labels;

    % Save mean zscore 
    ganalysis{pp,ii}.meanzsc=mean(data.zs,3) ;

end

%make sure that all Time vectors have the same length
%Get the lower length of alla ganalysis.Time vectors
tmp = [ganalysis{:}] ; 
minL = min(cell2mat(cellfun(@(x)(length(x)), {tmp.Time}, 'UniformOutput', false))) ; 

% Remove extra Time point if needed
for pp=1:length(ganalysis)
%     for ii=1:length(sFiles) 
        ganalysis{pp,1}.Time=ganalysis{pp,ii}.Time(1:minL);
        ganalysis{pp,1}.tvals=ganalysis{pp,ii}.tvals(:,1:minL);
        ganalysis{pp,1}.pvals=ganalysis{pp,ii}.pvals(:,1:minL);
%     end
end


%Close waitbar
close(hwait)
