function [ganalysis] = s5_group_data_v1(sFiles,mia_table,OPTIONS)
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
% Copyright (C) 2016-2018 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)

% analysis directory (must contain .eeg .vhdr and .vmrk)

ganalysis=[];

[subj,b,~]=unique(mia_table(:,1)); %subj = list of all the selected subject
%b=line number of the first file
%coresponding to another subject
count=cat(1,b,length(sFiles)+1);  %count will be used to count number of selected files owning by the same subject and put it in the rigth line
%  (should not be usefull because of the test existing on lines 33 to 43 )

%waitbar initialization
hwait=waitbar(0,sprintf('%s %d/%d',strrep(char(subj(1)),'_',' '),1,length(subj)),'Name','Group analysis');

% Make the waitbar stay on top
set(hwait,'WindowStyle','modal')

for pp=1:length(subj)
    
    waitbar(pp/length(subj),hwait,sprintf('%s %d/%d',strrep(char(subj(pp)),'_',' '),pp,length(subj)))
    
    %column counter
    n=1;
    %Find the right data and stats name
    for ii=count(pp):(count(pp+1)-1)
        [PATHSTR,filename,EXT1] = fileparts(sFiles{ii});
        [PATHSTRtmp,patient,EXT] = fileparts(char(PATHSTR));
        
        
        %check for LFP signal if zs has not been calculated yet
        if isempty(strfind(filename,'data'))
            sFile=fullfile(PATHSTR,strcat(filename,'_data.mat'));
            if ~exist(sFile,'file')
                % compute zs for LFP signal and give back the new sFile
                % name corresponding to zs
                zsOPTIONS.overwrite='yes';
                sFile=mia_s3_zscore(sFiles(ii),zsOPTIONS);
            end
        else
            sFile=sFiles{ii};
        end
       
        %load the sFile data (F, Time, ZS, freqb)
        data=load(char(sFile));
        
        %keep the intressing datas in ganalysis
        ganalysis{pp,n}.fname=sFile;
        ganalysis{pp,n}.subj=patient;
        ganalysis{pp,n}.df=size(data.F,3);
        ganalysis{pp,n}.freqb=data.freqb;
        ganalysis{pp,n}.edges=[OPTIONS.low_twin, OPTIONS.up_twin];
        
        
        %create the 'stats' file name that should exist
        statfile = strrep(sFile,'_data','_stats');
        %computation to find or calculate threshdur
        if exist(char(statfile))%check if stat file already exist
            stat_data=load(statfile); %load it
            if isfield(stat_data,'stats') %check for the existance of the field 'stats' (case of 'stat' file not computed by marspower)
                for ss=1:length(stat_data.stats)
                    if ~isempty(stat_data.stats(ss).fname)&& ~isempty(stat_data.stats(ss).pthresh)&& ~isempty(stat_data.stats(ss).nboot)&& ~isempty(stat_data.stats(ss).threshdur)
                        if stat_data.stats(ss).nboot==OPTIONS.nboot && stat_data.stats(ss).pthresh==OPTIONS.alpha %check for the right combination of nboot and alpha
                            ganalysis{pp,n}.threshdur=stat_data.stats(ss).threshdur;
                            ganalysis{pp,n}.threshp=OPTION.alpha;
                            
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
                            ganalysis{pp,n}.threshdur=stat_data.stats(ss).threshdur; %keep the new calculated threshdur
                        end
                    end
                end
            end
            
        else  %cas of inexistent statfile
            statfile=mia_s5_compute_stats(cellstr(sFile),OPTIONS);
            stat_data=load(char(statfile));
            ganalysis{pp,n}.threshdur=stat_data.stats.threshdur;
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
        ganalysis{pp,n}.Time=data.Time(find(data.Time));
        ganalysis{pp,n}.tvals=tvals(:,find(data.Time));
        ganalysis{pp,n}.pvals=pvals(:,find(data.Time));
        
        % Place the p(for left contacts) at the right place on contact name
        data.labels=strrep(data.labels,'''','p');
        ganalysis{pp,n}.labels=data.labels;
        
        % Save mean zscore 
        ganalysis{pp,n}.meanzsc=mean(data.zs,3) ;
       
        n=n+1;
    end
 
end

%make sure that all Time vectors have the same length
%Get the lower length of alla ganalysis.Time vectors
MinL=length(ganalysis{1,1}.Time);
for pp=1:length(subj)
    n=1;
    for ii=count(pp):(count(pp+1)-1) %keep the same way to count columns to avoid possible empty cell array in ganalysis => ii counts for the number of ganalysis of each column without empty cell
        MinL=min(length(ganalysis{pp,n}.Time),MinL);
        n=n+1;
    end
end
%remove extra Time point if needed
for pp=1:length(ganalysis)
    n=1;
    for ii=count(pp):(count(pp+1)-1)
        Time=ganalysis{pp,n}.Time;
        tvals=ganalysis{pp,n}.tvals;
        pvals=ganalysis{pp,n}.pvals;
        Time=Time(1:MinL);
        tvals=tvals(:,1:MinL);
        pvals=pvals(:,1:MinL);
        ganalysis{pp,n}.Time=Time;
        ganalysis{pp,n}.tvals=tvals;
        ganalysis{pp,n}.pvals=pvals;
        n=n+1;
    end
end

waitbar(length(subj)/length(subj),hwait,sprintf('%s %d/%d',strrep(char(length(subj)),'_',' '),length(subj),length(subj)))

%Close waitbar
close(hwait)











