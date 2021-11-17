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
% ***********************************************************************
%
% Script to contrast two conditions :
%
% Pre-requirement : 
% 1) Both conditions should be processed separately within MIA
% 2) The result is two MIA database : MIA_db_Condition1, MIA_db_Condition2
%
% The script then :
% 1) Loads the data from the two conditons 
% 2) Compute the permutations 
% 3) Compute the statistics within each roi 
% 4) Compute FDR correction and display results 
% -----------------------------------------------------------------------

%% Load MIA's working directories per condition
cOPTIONS.indir1 = '/Path/of/MIA_db_Condition1';
cOPTIONS.indir2 = '/Path/of/MIA_db_Condition2';
cOPTIONS.mtg = 'bipolar';

% rOPTIONS.mtg = 'monopolar';
cOPTIONS.freq = cell2mat(grpOPTIONS.freqid(getOPTIONS.freq));
cOPTIONS.win_noedges = [-401,1601] ;

%% Prepare the signals for statistics (Get the trials signals when needed)
tic ; [cfroi] = mia_browse_csignal(froi,cOPTIONS); toc ; 

%% Compute statistical thrersholds (permutation based) : 
% Warning :  This may take a while 
pOPTIONS.win_noedges = [-401,1601] ;
pOPTIONS.threshp = 0.05;
pOPTIONS.smoth = 20 ;
pOPTIONS.nperm =1000 ;
tic ; [proi] = mia_stats_permutations_rois(cfroi,pOPTIONS);toc

%% Save the results for future display (recommended because previous steps are so long to process)
save('proi','proi'); 

%% Compute the statistics
dOPTIONS.clr = jet(numel(unique(m_table_as(:,1)))) ; 
dOPTIONS.win_noedges = [-401,1601] ;
dOPTIONS.threshp = 0.05;
dOPTIONS.smoth = 20 ; 
[lroi] = mia_rois2pvalues(proi,dOPTIONS);

%% Compute FDR and display thresholded tvalues (uncorrected, corrected duration, corrected sumtval)
dOPTIONS.clr = jet(numel(unique(m_table_as(:,1)))) ; 
dOPTIONS.win_noedges = [-401,1601] ;
dOPTIONS.threshp = 0.05;
dOPTIONS.smoth = 20 ; 
dOPTIONS.tw = [0,1601] ;
mia_display_rois_conditions_fdr(lroi,dOPTIONS);