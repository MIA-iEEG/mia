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
% Function that loads signals from two conditions (indir1, indir2) for
% specific regions and return the signals (sig1, sig2) into the same roi 
% structure
% -----------------------------------------------------------------------
function [roi] = mia_browse_csignal(roi,OPTIONS) 

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

sig1 = cell(size(m_tab,1),1);
sig2 = cell(size(m_tab,1),1);

for kk=1:length(un)

    tmp1 = dir(fullfile(OPTIONS.indir1,un{kk},strcat('*',OPTIONS.mtg,'*data_',num2str(OPTIONS.freq),'_','*.mat')));
    tmp2 = dir(fullfile(OPTIONS.indir2,un{kk},strcat('*',OPTIONS.mtg,'*data_',num2str(OPTIONS.freq),'_','*.mat')));

    fname1 = fullfile(OPTIONS.indir1,un{kk},tmp1.name);
    fname2 = fullfile(OPTIONS.indir2,un{kk},tmp2.name);

    dat1 = load(fname1);      
    dat2 = load(fname2);      
    
    % Get channel to load
    idx_elec_tab = ismember(m_tab(:,1),un{kk}) ;
    lab_toproc = m_tab(idx_elec_tab,2);
    idx_elec_dat = ismember(dat1.labels,lab_toproc);  % same labels dat1 or dat2 
    id = find(idx_elec_dat);
  
    [LIA,LOCB] = ismember(lab_toproc',dat1.labels(idx_elec_dat)) ;
    
    % Format output
    sig1(idx_elec_tab) =  num2cell(dat1.zs(id(LOCB),(dat1.Time*1000>OPTIONS.win_noedges(1))&(dat1.Time*1000<OPTIONS.win_noedges(2)),:),[2,3]);
    sig2(idx_elec_tab) =  num2cell(dat2.zs(id(LOCB),(dat2.Time*1000>OPTIONS.win_noedges(1))&(dat2.Time*1000<OPTIONS.win_noedges(2)),:),[2,3]);
%     idpt(idx_elec_tab) = repmat({repmat(kk,sum(isGoodinRt),1)},sum(idx_elec_tab),1);

end

% Reinject signal into structure 
for ii=1:length(roi)

    idx = [m_tab{:,3}] ==ii;
    roi{ii}.sig1 = sig1(idx);
    roi{ii}.sig2 = sig2(idx);
    
end
