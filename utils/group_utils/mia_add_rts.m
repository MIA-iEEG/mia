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

    pt_name = [] ; 
    
    % Get current ROI
    croi = roi{ii};
    
    % Quick drity fix to ghandle underscore in the patient name. 
    % Get patients in the ROI 
    for pp=1:length(croi.labels)
        chan = cell2mat(croi.labels(pp)) ; 
        k=strfind(chan,'_');
        pt_name{pp} =  chan(1:k(2)-1) ;
    end
    
    % Alternatively (no underscore in patient ID) use :
%     pt_name = cellfun(@(x) x(1:strfind(x,'_')-1), croi.labels, 'UniformOutput',false) ;
   
    % Just add rts if there were existing ones (WARNING existing ones will
    % be overwritten) 
    if ~isfield(croi,'rts') ; croi.rts = cell(size(pt_name,1),1); end 
    
    % Gets all unique occurence of patients in this roi 
    [lia1,locb1] = ismember(pt_name,pt) ; 
    
    idx = locb1(locb1~=0) ; 
    
    % Gets Rts for this ROIs 
    croi.rts(locb1~=0) =  {tmp(idx).rt} ; %num2cell(cell2mat([tmp(idx).rt]),1)';
            
    % Reinject current roi into struct array
    roi{ii} = croi ;
end
