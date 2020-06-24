function [struct_table, status, message] = read_loc_table(filename, OPTIONS)
% -------------------------------------------------------------------------
% DESCRIPTION
%   Reads localization table
%
% Inputs :
%         filename : Excel file containg the table formatted as follow :

%
% Output:   status  : integrity of table (-1 if doublons exist ; 1
% otherwise) 
%           message : string output message containing doublons (per patient) 
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
% Copyright (C) 2016-2020 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)
% The excel file should contain four columns. Each column should start with
% one of the keywords "Patient", "Electrodes", "Lateralization" and "Region"

%Load excel file
[~, text_data, all_data] = xlsread(filename) ;
text_data = all_data ;


%%
% If there is not exactly 4 col
if size(text_data,2)~=4
    struct_table = [] ; status = -1; message =  'Wrong format (table should contain 4 columns)' ;
%     errordlg('Wrong format (table should contain 4 columns)') ;
   return;
end

status = 1 ;

message = '' ; % init
struct_table = [];

% Get the header line and remove spaces
hdr = deblank(text_data(1,:)) ;
idpt = find(strcmpi(hdr,'Patient'));

% idelec can either be Electrode or Contact
idelec = find(strcmpi(hdr,'Electrode'));
if isempty(idelec) ; idelec = find(strcmpi(hdr,'Contact')); end

idlat = find(strcmpi(hdr,'Lateralization'));
idroi = find(strcmpi(hdr,'Region'));

if isempty(idpt)|| isempty(idelec) || isempty(idlat)|| isempty(idroi)
    status =-1; 
    message = 'No header was found in file';
else

% Removes blanks : if no patient ID removes the whole line
text_data(strcmp(text_data(:,1),''),:) = [] ;
text_data(strcmp(text_data(:,2),''),:) = [] ;
text_data(strcmp(text_data(:,3),''),:) = [] ;
text_data(strcmp(text_data(:,4),''),:) = [] ;

% Removes the header line
text_data = text_data(2:end,:);

ptname = text_data(:,idpt) ; 
ptname = strrep(ptname,'''',''); ptname = strrep(ptname,',','');  ptname = strrep(ptname,'"','');

% Get unique occurence of patients 
u_pt = unique(ptname) ;
    
% For each patient 
for iPt=1:length(u_pt) 
   
    % Get indice of this patient in the table
     idx =  strcmp(u_pt{iPt},ptname);
    
    % Get pt name
    pt = u_pt{iPt};
    lat = text_data(idx,idlat);
    roi = text_data(idx,idroi);
    elec = text_data(idx,idelec);
        
    % Converts numeric region name (if any)
    idx = cell2mat(cellfun(@isnumeric, roi,'UniformOutput', false));
    roi(idx) = cellstr(num2str(cell2mat(roi(idx)))) ;
    
    % Clean text (no ' or coma or anything else but char and digits)
    lat = strrep(lat,'''',''); lat = strrep(lat,',','');  lat = strrep(lat,'"','');
    roi = strrep(roi,'''',''); roi = strrep(roi,',',''); roi = strrep(roi,'"','');
    elec = strrep(elec,'''',''); elec = strrep(elec,',',''); elec = strrep(elec,'"','');
    
    % Create contacts labels (electrode+number+laterality) 
    contacts = strcat(elec,'(',lat,')');
   
    [a,ia1,ic1] = unique(contacts) ; 
     n = histc(ic1,unique(ic1));
    
    % Doublons were found
    doublons = a(find(n>=2)) ; 
    if ~isempty(doublons)
        listd = sprintf('%s, ',doublons{:});
        message = sprintf('%s\n %s : %s',message,pt,listd);
        status = -1;
    end
    
    struct_table{iPt}.pt = pt; 
    struct_table{iPt}.lat = lat ; 
    struct_table{iPt}.elec = elec; 
    struct_table{iPt}.roi= roi; 
    
end

end
