function [struct_table] = read_rt_table(filename, OPTIONS)
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
% Copyright (C) 2016-2018 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)

%Load excel file
[~, text_data, all_data] = xlsread(filename) ;
% %%
% % If there is not exactly 4 col
% if size(text_data,2)~=4
%     errordlg('Wrong format (table should contain 4 columns)') ;
%    return;
% end

struct_table = [];

% Get the header line and remove spaces
hdr = deblank(text_data(1,:)) ;
idpt = find(strcmpi(hdr,'Patient'));
idor = find(strcmpi(hdr,'Order'));
idrt = find(strcmpi(hdr,'RT'));
idat = find(strcmpi(hdr,'Indata'));

if isempty(idpt)|| isempty(idor) || isempty(idrt)|| isempty(idat)
    status =-1; 
    message = 'No header was found in file';
else

% Removes the header line
all_data= all_data(2:end,:);

% Removes blanks : if no patient ID removes the whole line
all_data(strcmp(all_data(:,idpt),''),:) = [] ;
all_data(isnan(cell2mat(all_data(:,idor))),:) = [] ;
all_data(isnan(cell2mat(all_data(:,idrt))),:) = [] ;
all_data(isnan(cell2mat(all_data(:,idat))),:) = [] ;


ptname = all_data(:,idpt) ; 
ptname = strrep(ptname,'''',''); ptname = strrep(ptname,',','');  ptname = strrep(ptname,'"','');

% Get unique occurence of patients 
u_pt = unique(ptname) ;
    
% For each patient 
for iPt=1:length(u_pt)  
   
    % Get indice of this patient in the table
     idx =  strcmp(u_pt{iPt},ptname);
    
    % Get pt name
    pt = u_pt{iPt};
    rt = all_data(idx,idrt);
    ord = all_data(idx,idor);
    dat= all_data(idx,idat);

    struct_table{iPt}.pt = pt; 
    struct_table{iPt}.rt = rt ; 
    struct_table{iPt}.order = ord; 
    struct_table{iPt}.indat= dat; 
    
end

end
