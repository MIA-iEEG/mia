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
function [roi] = mia_add_rts(roi,rt_table, OPTIONS) 

% Get RT table
tmp = [rt_table{:}];
pt = {tmp.pt}' ; 

%  Loop through regions
for ii=1:length(roi) 
   
    % Get current ROI
    croi = roi{ii};
    
    % Just add rts if there were existing ones (WARNING existing ones will
    % be overwritten) 
    if ~isfield(croi,'rts') ; croi.rts = cell(length(croi.labels),1); end 

    % Find patient name in electrodes labels (if RTs exist for this
    % patient, otherwise '')
    for pp=1:length(croi.labels)
        if sum(cellfun(@(x) contains(croi.labels(pp),x),pt))~=0
            pt_name{pp} = pt{cellfun(@(x) contains(croi.labels(pp),x),pt)} ; 
        else 
            pt_name{pp} = '';
        end
        
    end
      
    % Gets all unique occurence of patients in this roi 
    [lia1,locb1] = ismember(pt_name,pt) ; 
    
    idx = locb1(locb1~=0) ; 
    
    % Gets Rts for this ROIs 
    croi.rts(locb1~=0) =  {tmp(idx).rt} ; %num2cell(cell2mat([tmp(idx).rt]),1)';
            
    % Reinject current roi into struct array
    roi{ii} = croi ;
end
