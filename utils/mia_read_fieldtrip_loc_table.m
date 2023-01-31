function [struct_table, status, message] = mia_read_fieldtrip_loc_table(filename, OPTIONS)
% -------------------------------------------------------------------------
% DESCRIPTION
%   Reads a Brainstorm iEEG atlas table (.tsv)
%
% Inputs :
%         filename : Name of the iEEG atlas table (filename.tsv)

%
% Output:   status  : integrity of table (-1 if doublons exist ; 1
% otherwise) 
%           message : string output message containing doublons
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
% Copyright (C) 2016-2022 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)

%Load excel file
[~, text_data, all_data] = xlsread(filename) ;
text_data = all_data ;

status = 1 ;

message = '' ; % init
struct_table = [];

% Get the header line and remove spaces
hdr = deblank(text_data(1,:)) ;
idpt = find(strcmpi(hdr,'Patient'));


% Init output variables 
message = '' ; 
struct_table = [];
status = 1 ; 

% Get the header line and remove spaces
hdr = deblank(text_data(1,:)) ;

% idelec can either be Electrode or Contact
idelec = find(strcmpi(hdr,'Electrode'));

not_atlas = {'Electrode','Coordinates','Discard','Epileptic', 'Out of Brain', 'Notes','Loc Meeting'} ; 

atlases = hdr(~ismember(hdr,not_atlas)) ; 

% Something is wrong with the format (missing column or missing header)
if isempty(idelec) 
    status =-1; 
    message = ' Table should contain a column "Electrode"';
else
  
    % Converts numeric to characters
    idx_numeric=cellfun(@isnumeric,text_data,'uni',false); idx_numeric = cellfun(@any,idx_numeric) ; 
    text_data(idx_numeric) = cellfun(@num2str, text_data(idx_numeric),'uni',false) ;

    % Removes the header line
    text_data = text_data(2:end,:);

    [ptname,atlas]=mia_inputdialog(OPTIONS.patients,atlases) ; 

    if isempty(ptname)||isempty(atlas) ; status =2 ; return ; end 
   
    % Find columns containg regions labels
    idroi = find(strcmpi(hdr,atlas));
    roi = text_data(:,idroi);
    
    % Get laterality from coordinates (x>0 = Rigth ; x<0 = Left)
    idcoord = find(strcmpi(hdr,'Coordinates'));
    coord  = text_data(:,idcoord) ;
    idx_left = find(strcmp(cellfun(@(x) x(1), coord, 'UniformOutput',false),'-')) ; 
    idx_right = find(~strcmp(cellfun(@(x) x(1), coord, 'UniformOutput',false),'-')) ; 
  
    lat(idx_left) = {'L'}; lat(idx_right) = {'R'}; 
    
    % Get electrode names
    idelec = find(strcmpi(hdr,'Electrode'));
    elec = text_data(:,idelec);
    
    % Removes any anoying character
    elec = strrep(elec,'''',''); elec = strrep(elec,',',''); elec = strrep(elec,'"','');

    % There is only patient in a file that is why we use {1} (otherwise
    % several patient table can be stored in struct_table
    struct_table{1}.pt = ptname; 
    struct_table{1}.lat = lat ; 
    struct_table{1}.elec = elec; 
    struct_table{1}.roi= roi; 
    struct_table{1}.atlas= atlas; 

end
end

function [opt1,opt2]=mia_inputdialog(opt_list1,opt_list2)

% Create dialog box with two scroll down menus 
hfig=figure('CloseRequestFcn',@close_req_fun,'menu','none','Units','Normalized', ...
    'Position', [.5, .5, .1, .1],...
    'Name','Select patient and atlas',...
    'menu','none',...
    'NumberTitle','off');

opt1= [] ; opt2  =[] ;

dropdown1=uicontrol('Style', 'popupmenu', 'String', opt_list1, ...
    'Fontsize',12,...
    'Parent',hfig,'Units','Normalized', ...
    'Position',  [.1, .65, .8, .15]);
dropdown2=uicontrol('Style', 'popupmenu', 'String', opt_list2, ...
    'Fontsize',12,...
    'Parent',hfig,'Units','Normalized', ...
    'Position', [.1, .4, .8, .15]);
uicontrol('Style', 'pushbutton', 'String', 'OK', ...
    'Parent',hfig,'Units','Normalized', ...
    'Fontsize',12,...
    'Position', [.1 .1 .35 .2],...
    'Callback','close(gcbf)');
cancel=uicontrol('Style', 'pushbutton', 'String', 'cancel', ...
    'Fontsize',12,...
    'Parent',hfig,'Units','Normalized', ...
    'Position', [.55 .1 .35 .2],...
    'Tag','0','Callback',@cancelfun);

% Wait for figure being closed (with OK button or window close)
uiwait(hfig)

% Figure is now closing
if strcmp(cancel.Tag,'0')%not canceled, get actual inputs
    opt1=opt_list1{dropdown1.Value};
    opt2=opt_list2{dropdown2.Value};
end

% Actually close the figure
delete(hfig)

end
function cancelfun(h,~)
set(h,'Tag','1')
uiresume
end
function close_req_fun(~,~)
uiresume
end
