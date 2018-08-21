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
function [rasters,idpt,m_tab] = get_rasters(roi,OPTIONS) 

ct=1;
%  Loop through regions
for ii=1:length(roi) 

    croi = roi{ii};
    
    k=strfind(croi.labels,'_');
    
    for pp=1:length(k)
    
        str = cell2mat(croi.labels(pp)) ;
        
        % Column for pt ID
        if strcmpi(OPTIONS.mtg,'monopolar') ; tt = 0 ; else tt=1 ; end
        
        m_tab{ct,1} = str(1:k{pp}(end-tt)-1);
        
        % Column for electrode
        m_tab{ct,2} = str(k{pp}(end-tt)+1:end);
        
        % Add column for indexing the roi
        m_tab{ct,3} = ii;
        ct=ct+1;
        
    end
end

% Gets all unique occurence of patients
[un, ia, ic] = unique(m_tab(:,1),'stable');

rasters = cell(size(m_tab,1),1);
rts = cell(size(m_tab,1),1);

for kk=1:length(un)

    tmp = dir(fullfile(OPTIONS.maindir,un{kk},strcat('*',OPTIONS.mtg,'*data*',num2str(OPTIONS.freq),'*.mat')));

    fname = fullfile(OPTIONS.maindir,un{kk},tmp.name);

    fileg = dir(fullfile(OPTIONS.maindir,un{kk},'*signal_LFP*.mat'));
          
    dat = load(fname);      
    
    % Get channel to display
    idx_elec_tab = ismember(m_tab(:,1),un{kk}) ;
    lab_toproc = m_tab(idx_elec_tab,2);
    idx_elec_dat = ismember(dat.labels,lab_toproc);    
    
%     [C,IA,IC] = unique(lab_toproc,'stable');
    
    id = find(idx_elec_dat);
  
    [LIA,LOCB] = ismember(lab_toproc',dat.labels(idx_elec_dat)) ;
    
    % Format output
    rasters(idx_elec_tab) =  num2cell(dat.zs(id(LOCB),(dat.Time>OPTIONS.win_noedges(1))&(dat.Time<OPTIONS.win_noedges(2)),:),[2,3]);
    idpt(idx_elec_tab) = repmat({repmat(kk,size(dat.zs,3),1)},sum(idx_elec_tab),1);

end

