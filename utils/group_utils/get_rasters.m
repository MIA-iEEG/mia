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
function [roi] = get_rasters(roi,OPTIONS) 

ct=1;
%  Loop through regions
for ii=1:length(roi) 

    croi = roi{ii};
    
    k=strfind(croi.labels,'_');
    
    for pp=1:length(k)
    
        str = cell2mat(croi.labels(pp)) ;
        
        % Select rigth index of the '_' if monopolar and/or Flipped (_FLP)
        if strcmpi(OPTIONS.mtg,'monopolar')&&length(k{pp})<=1 ; tt = 0 ; else tt=1 ; end
        
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
idpt = cell(size(m_tab,1),1);

% Create progress bar for patients processing
hwait_pt = waitbar(0,'','Units','Normalized','Name','Loading data...');
% Make the waitbar stay on top
set(hwait_pt,'WindowStyle','modal')

% For each patient load single trial data 
for kk=1:length(un)

    % Update progress bar 
    waitbar(kk/length(un),hwait_pt,sprintf('%s %s','Loading rasters',un{kk})) ;

%     tmp = dir(fullfile(OPTIONS.maindir,un{kk},strcat('*',OPTIONS.mtg,'*data*',num2str(OPTIONS.freq),'.mat')));
    tmp = dir(fullfile(OPTIONS.maindir,un{kk},strcat('*',OPTIONS.mtg,'*data*','.mat')));
    
    % A montage file was found
    if size(tmp,1)~=1
        tmp = tmp(cell2mat(cellfun(@isempty, strfind({tmp.name},'montage'), 'UniformOutput',false)));
    end
    
    fname = fullfile(OPTIONS.maindir,un{kk},tmp.name);

    fileg = dir(fullfile(OPTIONS.maindir,un{kk},'*signal_LFP*.mat'));
          
    dat = load(fname);      
   
    % Replace all ' by p
    dat.labels = strrep(dat.labels,'''','p') ; 
    
    % Get channel to display
    idx_elec_tab = ismember(m_tab(:,1),un{kk}) ;
    lab_toproc = m_tab(idx_elec_tab,2);
    idx_elec_dat = ismember(dat.labels,lab_toproc);    
    
    % Define time vector once for all 
    if kk==1; vTime = (dat.Time>roi{1}.t(1))&(dat.Time<=roi{1}.t(end)) ; end
        
    id = find(idx_elec_dat);
  
    [LIA,LOCB] = ismember(lab_toproc',dat.labels(idx_elec_dat)) ;
    
    % TODO HERE : debug won't work if there are some flipped channels...
    % (look at 2nd colum in m_tab : contact labels including _FLP wont get read)
    % Format output
    rasters(idx_elec_tab) =  num2cell(dat.zs(id(LOCB),vTime,:),[2,3]);
    idpt(idx_elec_tab) = repmat({repmat(kk,size(dat.zs,3),1)},sum(idx_elec_tab),1);
    % TODOS : ASD add indat information (isGood -> contained in data) 
    

end
  
% Reinject the signals into the roi structure 
for ii=1:length(roi) 

    idx = [m_tab{:,3}] == ii ;   
    roi{ii}.F =rasters(idx);
    roi{ii}.idptsignals= idpt(idx);

end

% Load is done : close progress bar
delete(hwait_pt) ;
