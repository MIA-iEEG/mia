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
%
% This function browse recursivvely in working directory to fill
% out a table containgin all files (processes) that have been produced
% 2017/1/17 : ASD add options.allow_signflip
function [rois,mia_table,edges,datafiles]=create_struct_of_rois(m_table_as,ganalysis, OPTIONS)

isGood =~(strcmp(m_table_as(:,5),'out')|strcmp(m_table_as(:,5),'lesion')|strcmp(m_table_as(:,5),'les'));
% [m_table_effect, s, smask, all_labels] = mia_get_table_effect_clc(m_table_as(isGood,:), ganalysis);
[m_table_effect, s, smask, all_labels] = mia_get_table_effect_clc(m_table_as(isGood,:), ganalysis);

%Define OPTIONS for group anlaysis and display
% getOPTIONS.freq = 1; % Morlet or Hilbert what freqid? 
getOPTIONS.freq = 1; % Morlet or Hilbert what freqid? 
getOPTIONS.nPt= 1 ;% min numb of pt by roi
getOPTIONS.signifmode =1 ; % mode to select significant activity (if 0 no constrain)
getOPTIONS.signmode = 'signed' ;

% getOPTIONS.signmode = 'abs';

gan = [ganalysis{1,:}] ;
 
% if flip sign checked
if OPTIONS.allow_flipsign 
    getOPTIONS.flip_thresh = OPTIONS.flip_thresh ; 
end

% Get all rois that has significant activity for this frequency band
rois = mia_get_roi(m_table_effect,ganalysis{1,1}.Time, s, smask, all_labels,{gan.freqb},getOPTIONS);

%remove blanc
r = [rois{:}];
rois(strcmp({r.name},'blancL')|strcmp({r.name},'blancR'))=[];

id_name=1;
id_onset=2;
id_corrPt=3;
id_corrChan=4;
id_numPT=5;
id_ID=6;

mia_table = [] ; 

% Fill table with rois 
for ii=1:length(rois)
    
    mia_table{ii,id_name}=rois{ii}.name;
    mia_table{ii,id_onset}=num2str(rois{ii}.onset);
    mia_table{ii,id_corrPt}=num2str(rois{ii}.corrPt);
    mia_table{ii,id_corrChan}=num2str(rois{ii}.corrChan);
    mia_table{ii,id_numPT}=num2str(length(rois{ii}.namePt));
    mia_table{ii,id_ID}=num2str(ii);
    
end

%get back time window information 
edges=ganalysis{1}.edges;

%get back list of data files for this study 
n=1;
for ii=1:size(ganalysis,1)
    for jj=1:size(ganalysis,2)
        if ~isempty(ganalysis{ii,jj})
         datafiles{n,1}=char(ganalysis{ii,jj}.fname);
         n=n+1;
        end
    end
end











