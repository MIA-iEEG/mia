function [struct_table, status, message] = mia_read_loc_tsv_table(filename, OPTIONS)
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

% Init output variables 
message = '' ; 
struct_table = [];
status = 1 ; 

% Read the tsv file
T = readtable(filename,'FileType','text','ReadVariableNames',true) ;
% This throw a warning message due to semi column present in some headers 
% (e.g. cortex_148917V:Schaefer_100_17net) Matlab replaces with _
% 

% Get table column headers
atlases = T.Properties.VariableNames ;

% Excludes _prob and coordinates
atlases = atlases(~contains(atlases,'_prob')) ;
atlases = atlases(~ismember(atlases,'Channel')&...
                    ~ismember(atlases,'SCS')&...
                    ~ismember(atlases,'MNI')&...
                    ~ismember(atlases,'World')) ; 

[ptname,atlas]=mia_inputdialog(OPTIONS.patients,atlases) ; 

if isempty(ptname)||isempty(atlas) ; status =2 ; return ; end 
    
idx = find(ismember(fieldnames(T),atlas)) ; 

% Get region labels
roi = T.(idx);

idx_left = contains(roi,' L') | contains(roi,'Left') ; 
idx_right = contains(roi,' R') | contains(roi,'Right') ; 

lat(idx_left) = {'L'}; lat(idx_right) = {'R'}; 

% Only keeps data for which we have a lateraltiy (exclude de facto N/A)
roi = roi(idx_left|idx_right);
lat = lat(idx_left|idx_right)';
elec = T.Channel(idx_left|idx_right);

% Removes any anoying character
elec = strrep(elec,'''',''); elec = strrep(elec,',',''); elec = strrep(elec,'"','');

% % Something is wrong with the format (missing column or missing header)
% if isempty(idpt)|| isempty(idelec) || isempty(idlat)|| isempty(idroi)
%     status =-1; 
%     message = ' Table should contain 4 columns : "Patient", "Contact (or Electrode)", "Lateralization","Region" ';
% else

struct_table{1}.pt = ptname; 
struct_table{1}.lat = lat ; 
struct_table{1}.elec = elec; 
struct_table{1}.roi= roi; 
struct_table{1}.atlas= atlas; 

% end
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