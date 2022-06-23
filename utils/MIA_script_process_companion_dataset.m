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
% Copyright (C) 2016-2022 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)
% ***********************************************************************
%
% Script to analyse intracranial dataset provided on Zenodo 
% (doi 10.5281/zenodo.4767855)
% 
% 1) Load the Zenodo companion dataset 
% 2) Create a MIA database 
% 3) Compute a time-frequency decomposition on all patients 
% 4) Compute non-parametric statistics on each patient 
% 5) Loads the companion anatomical labeling table
% 6) Creates a MIA study (which can be opened from the GUI : 
%       "Display" button top rigth corner under "Studies" panel) 
%
% -----------------------------------------------------------------------

%% BUILD MIA out of Zenodo files
%-----------------------------------------------------------------------------

%% CREATE directory structure
% Path to the directory in which you downloaded ALL Zenodo data files
zenodo_import_dir = '/Users/anne-sophiedubarry/Documents/projects/in_progress/MIA/data/Zenodo_companion_data';

% Read data files from Zenodo local directory
d = dir(zenodo_import_dir); 
idx = ~cellfun(@isempty, strfind({d.name},'_signal_LFP.mat')) ; 
datafiles= {d(idx).name} ; 
patients_dir = strrep(datafiles, '_signal_LFP.mat','') ; 

% Create a MIA_db directory in zenodo' folder
mia_db = fullfile(zenodo_import_dir,'MIA_db') ; 
if ~exist(mia_db, 'file') ; mkdir(mia_db) ; end

% Move all patient data files into one folder patient directory 
for pp=1:length(patients_dir)
    pt_dir = fullfile(mia_db,patients_dir{pp}) ; 
    if ~exist(pt_dir, 'dir') 
        mkdir(pt_dir);
        copyfile(fullfile(zenodo_import_dir,datafiles{pp}),pt_dir);
    end
end

%% Decompose Time-Frequency 
%-----------------------------------------------------------------------------
%Options = Select descomposition mode, frequency band and baseline
OPTIONS.modetf= 'Morlet'; 
OPTIONS.mtg = 'BIPOLAR';
OPTIONS.outdir = mia_db ; 
OPTIONS.freqs= [80:10:150]; % lower freq, step, higher freq : [lfreq:step:ufreq]
OPTIONS.ncycles = 7 ; 
% OPTIONS.freqs= logspace(log10(80), log10(150),8) ; 
OPTIONS.subjects = patients_dir ; 
OPTIONS.removeEvoked = 0; % Do not remove evoked response prior to TF decomposition
OPTIONS.zbaseline = [ -0.9721, -0.0289] ; % Define baseline for nomalisation
OPTIONS.overwrite = 'No';

% Compute frequecy analysis
sFiles = mia_s4_compute_tf_bandwise(OPTIONS);

% Launch MIA (use Display button to check results)
mia(mia_db);

%% Compute statistics
%-----------------------------------------------------------------------------
groupOPTIONS.nboot=3;
groupOPTIONS.alpha=0.001;
groupOPTIONS.twin=[-0.5,1.5];
groupOPTIONS.outdir=mia_db;
groupOPTIONS.baseline=[ -0.8, -0.05] ;
groupOPTIONS.subjects = patients_dir ; 
    
% Check if GA_Results exists, if not create it
[upath,~,~]=fileparts(mia_db);
outdir=fullfile(upath,'GA_Results');
if ~exist(outdir,'dir') ; mkdir(outdir); end

% Set name of the study 
study_name = 'HGA_linearStep' ; 

% Check if this specific study exist
if exist(fullfile(outdir,study_name)) 
    error('Study exist, change study name to create a new one') ; 
else
    % Create study dir
    mkdir(char(fullfile(outdir,study_name)));
end

% Get one structure for group analysis (ganalysis)
ganalysis = mia_s5_group_data_v1(sFiles,groupOPTIONS);

% Create a study filename
[~,b] = fileparts(ganalysis{1}.fname) ; 
[mtg remain] = strtok(b(length(patients_dir{pp})+1:end),'_') ;
[method remain] = strtok(remain,'_') ;
[datatype remain] = strtok(remain,'_') ;
strfreq = strsplit(remain,'_');
freq = sprintf('%s-%s Hz (step %s)',strfreq{3},strfreq{4},strfreq{2});
                    
fname =char(fullfile(outdir,study_name,strcat(study_name,'_',num2str(length(patients_dir)),'_',OPTIONS.modetf,'_',OPTIONS.mtg,'_',freq)));
save(fname,'ganalysis');


%% Load labeling table
%-----------------------------------------------------------------------------
% Loads the localizations of all patients contacts 
[struct_table, status, message] = mia_read_loc_table(fullfile(zenodo_import_dir,'loc_table_MIA.xlsx'),groupOPTIONS) ;    

loctable_name = 'MyTable' ;
groupOPTIONS.maindir = mia_db ; 

% Map the contacts from loc table with the ones in the data 
[m_table_all, status, message] = mia_get_dataloc_table(struct_table,groupOPTIONS);
  
% Atlas filename 
fname = fullfile(outdir,strcat('m_table_',loctable_name,'.mat'));
s.table_fname = fullfile(zenodo_import_dir,'loc_table_MIA.xlsx') ;
    
% Save m_table fopr this atlas
save(fname,'-struct','s');

%% GET TABLE  OF EFFECTS 
%-----------------------------------------------------------------------------
% Filters out regions not interessting (out  blanc, etc.0)
isGood =~(strcmp(m_table_all(:,5),'out')|strcmp(m_table_all(:,5),'lesion'));
[m_table_effect, s, smask, all_labels] = mia_get_table_effect_clc(m_table_all(isGood,:), ganalysis);


%% GET ROIS No stats on correlation 
%-----------------------------------------------------------------------------
%Define OPTIONS for group anlaysis and display
getOPTIONS.freq = 1; % Morlet or Hilbert what freqid? 
getOPTIONS.nPt= 2 ;% min numb of pt by roi
getOPTIONS.signifmode =1 ; % mode to select significant activity (if 0 no constrain)
getOPTIONS.signmode = 'signed';
getOPTIONS.montage = 'bipolar';
% getOPTIONS.signmode = 'abs';

gan = [ganalysis{1,:}] ;
 
% Get all rois that has significant activity for this frequency band
rois = mia_get_roi(m_table_effect,ganalysis{1,1}.Time, s, smask, all_labels,{gan.freqb},getOPTIONS);


%% GET ROIS with stats on correlations *
%-----------------------------------------------------------------------------
%Define OPTIONS for group anlaysis and display
getOPTIONS.freq = 1; % Morlet or Hilbert what freqid? 
getOPTIONS.nPt= 1 ;% min numb of pt by roi
getOPTIONS.signifmode =0 ; % mode to select significant activity (if 0 no constrain)
getOPTIONS.signmode = 'signed';
% getOPTIONS.signmode = 'abs';
gan = [ganalysis{1,:}] ;
 
% Get all rois that has significant activity for this frequency band
rois = mia_get_roi_permute(m_table_effect,ganalysis{1,1}.Time, s, smask, all_labels,{gan.freqb},getOPTIONS);

%% DISPLAYS Statistics on correlations
%-----------------------------------------------------------------------------
dOPTIONS.clr = jet(numel(unique(m_table_all(:,1)))); % distinct colors for pt
dOPTIONS.win_noedges = [-0.5,1.4] ;
mia_display_roi_and_permcorr(rois,dOPTIONS);


%% DISPLAYS SUMMARY (MEAN ROIS)
%-----------------------------------------------------------------------------
pOPTIONS.thresh =3; % -1 for no color chronological organization
pOPTIONS.nsub =4; % number of subplot
pOPTIONS.title = '';
[labels_o, colorm] = mia_display_summary_roi(rois,pOPTIONS) ;

%% DISPLAYS ONE PATIENT DATA (MEAN ROIS)
%-----------------------------------------------------------------------------
pOPTIONS.thresh = 3; % -1 for no color chronological organization
pOPTIONS.nsub =2; % number of subplot
pOPTIONS.ptKey = 'PT_02';

[labels_o_pt, colorm_pt] =mia_display_summary_roi_pt(rois,pOPTIONS) ;
