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


% This function browse recursivvely in working directory to fill
% out a table containgin all files (processes) that have been produced
% 2017/1/17 : ASD add options.allow_signflip
% 2018/4/19 : ASD add output montage (for raster display) 
function [rois,mia_table,edges,datafiles, montage, freqb, wdir]=create_table_of_rois(m_table_as,ganaFile, OPTIONS)

% Loads study file 
gana=load(ganaFile,'ganalysis');
ganalysis=gana.ganalysis;

%% Get montage (from one of the filename : at this stage all have the same montage) 
[pathname,fname,~] = fileparts(ganalysis{1}.fname);
tmp = strsplit(fname(length(ganalysis{1}.subj)+2:end),'_') ;
montage = tmp{1} ;

% FIX old files : Get threshp wich was used to compute the threshdur
if ~isfield(ganalysis{1,1},'threshp')
    threshp = str2double(input('WARNING : you used an old version of MIA to generate the study, please indicate the alpha critical values that you used : alpha = ', 's'));
    while isnan(threshp) 
        threshp = str2double(input('alpha (must be a number): ', 's'));
    end
    
    % The following 4 lines add the field 'threshp' to each element of cell
    % array ganalysis (prevents from looping)
    tmp = cell2mat(ganalysis) ;
    T = num2cell(threshp*ones(1,length(ganalysis)));
    [tmp(1:length(ganalysis)).threshp] = T{:};    
    ganalysis = num2cell(tmp);
    
    save(ganaFile,'ganalysis');
end

% Create table of effects 
[m_table_effect, s, smask, all_labels] = get_table_effect_clc(m_table_as, ganalysis);

%Define OPTIONS for group anlaysis and display
getOPTIONS.freq = 1; % Morlet or Hilbert what freqid? 
getOPTIONS.nPt= 1 ;% min numb of pt by roi (set to 1 for GUI : user can filter the ROIs by visualizing the number of pts from the table)
getOPTIONS.signifmode =OPTIONS.signifmode ; % mode to select significant activity (if 0 no constrain)
getOPTIONS.signmode = 'signed' ; % signmode can be 'signed' or 'abs';
getOPTIONS.montage = montage ; 

if OPTIONS.allow_flipsign 
    getOPTIONS.flip_thresh = OPTIONS.flip_thresh ; 
end

% Get all rois that has significant activity for this frequency band
rois = get_roi(m_table_effect,ganalysis{1,1}.Time, s, smask, all_labels,{ganalysis{1}.freqb},getOPTIONS);

% mia_table column index 
id_name=1;
% id_onset=2;
id_corrPt=2;
id_corrChan=3;
id_numPT=4;
id_numCT=5;
id_ID=6;

mia_table = [] ; 
cPt = [] ; 
cChan = [] ; 

%% Fill in table with rois
for ii=1:length(rois)
    
    mia_table{ii,id_name}=rois{ii}.name;
%     mia_table{ii,id_onset}=num2str(rois{ii}.onset);
    if isnan(rois{ii}.corrPt); cPt = '-'; else ; cPt = num2str(rois{ii}.corrPt);end
    mia_table{ii,id_corrPt}=cPt ; 
    if isnan(rois{ii}.corrChan); cChan = '-'; else ; cChan = num2str(rois{ii}.corrChan);end
    mia_table{ii,id_corrChan}=cChan;
    mia_table{ii,id_numPT}=num2str(length(rois{ii}.namePt));
    mia_table{ii,id_numCT}=num2str(size(rois{ii}.Fmask,2));
    mia_table{ii,id_ID}=num2str(ii);
    
end

%% Get time window information 
edges=ganalysis{1}.edges;

%% Get list of data files for this study 
n=1;
for ii=1:size(ganalysis,1)
    for jj=1:size(ganalysis,2)
        if ~isempty(ganalysis{ii,jj})
         datafiles{n,1}=char(ganalysis{ii,jj}.fname);
         n=n+1;
        end
    end
end


%% Get freq 
freqb = ganalysis{1}.freqb;

%% Get working directory
wdir = fileparts(pathname);

