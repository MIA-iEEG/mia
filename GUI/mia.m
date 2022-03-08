function varargout = mia(varargin)
% MIA MATLAB code for mia.fig
%      MIA, by itself, creates a new MIA or raises the existing
%      singleton*.
% 
% This GUI was created with GUIDE
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
% Copyright (C) 2016-2021 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @mia_OpeningFcn, ...
    'gui_OutputFcn',  @mia_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before mia is made visible.
function mia_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for mia
handles.output = hObject;

% Get User home directory
dirname = mia_getuserdir ;

% Set initial default OPTIONS
handles.extOPTIONS.suffix = '';
handles.extOPTIONS.overwrite = 'no';
handles.extOPTIONS.SensorType = 'SEEG'; % {'MEG GRAD','SEEG','MEG MAG'};
handles.extOPTIONS.mtg= 'Bipolar';
handles.current_loctable = [];

% Move window to the center of the screen 
movegui(gcf,'center');

% Loads database (either from filepath in varagin, history file or promtp
% user
if  ~isempty(varargin)
    % BST plugin was used to call MIA with a database dir 
    directoryname = varargin{1} ; 
    
elseif exist(fullfile(dirname,'.mia_history.mat'),'file')

    % if history exist set the working directory 
    hist = load(fullfile(dirname,'.mia_history.mat')); 
    directoryname = hist.history.dirname ; 
    handles.history = hist.history ;
    set(handles.outdir,'String',hist.history.dirname);
    
else
        
     % Open a directory browser
    hwarn= warndlg('This is the first time you run MIA. Please pick an empty directory for MIA to store its database','Database');
    waitfor(hwarn);
    directoryname = uigetdir('Pick a DataBase Directory');
    
end

if directoryname~=0
    
    % Set the new working directory
    handles = set_workingdirectory(handles,directoryname);

    % Initialize the GUI
    handles = initialize_gui(hObject, handles, false);

    % Update handles structure
    guidata(hObject, handles);

else
    close(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = mia_OutputFcn(hObject, eventdata, handles)

%% SAVE HISTORY

% Get User home directory
dirname = mia_getuserdir ;

% If user started to use MIA
if isfield(handles,'history')
    
    % Save history file
    history = handles.history ; 
    save(fullfile(dirname,'.mia_history.mat'),'history');

    % Get default command line output from handles structure
    varargout{1} = handles.output;

end


% --------------------------------------------------------------------
function handles = initialize_gui(fig_handle, handles, isreset)

% initialize the whole interface 
handles = update_patientlist(handles) ;
handles = update_studies_list(handles) ;
handles = update_atlas_list(handles) ;
handles = create_data_table(handles) ; 
handles = update_loctablepath_text(handles);

% Create the table (list) of data on the central panel
function handles = create_data_table(handles)

% Get table that contains all files (names) found in the working directory
[handles.table.mia_table,handles.table.sFiles] = mia_create_table_workdir(handles.extOPTIONS.outdir, handles.current_loctable) ;

jtable = com.jidesoft.grid.SortableTable(handles.table.mia_table,{'Patient','Method','Montage','Freq. band','Fs','Remove Avg','Nb stats','Localized Contacts','ID'});

% Trick to hide to indexing column (in last column : jtable.getColumnCount - 1)
jtable.getColumnModel().getColumn(jtable.getColumnCount - 1).setMinWidth(0);
jtable.getColumnModel().getColumn(jtable.getColumnCount - 1).setMaxWidth(0);
jtable.getColumnModel().getColumn(jtable.getColumnCount - 1).setWidth(0);

% Create the java table
tableHeader = com.jidesoft.grid.AutoFilterTableHeader(jtable);
tableHeader.setAutoFilterEnabled(true)
tableHeader.setShowFilterName(true) 
tableHeader.setShowFilterIcon(true)
jtable.setTableHeader(tableHeader)

% % Save table components in handles
handles.table.jtable = jtable ;

handles.table.jScrollPane = javax.swing.JScrollPane(handles.table.jtable);

% Set the jtable uneditable 
sclass=java.lang.String('').getClass;
renderer=jtable.getDefaultRenderer(sclass);
editor=jtable.getDefaultEditor(sclass);
editor.getComponent.setEditable(0);
editor.setClickCountToStart(intmax);

% Set allcolumn uneditable
for iCol=0:size(handles.table.mia_table,2)-1 
    jtable.getColumnModel.getColumn(iCol).setCellRenderer(renderer)
    jtable.getColumnModel.getColumn(iCol).setCellEditor(editor);
end

% hjtable = handle(jtable,'CallbackProperties');
set(jtable,'MouseClickedCallback',{@edit_stats,handles});
set(jtable,'KeyReleasedCallback',{@edit_stats,handles});

% Copy the table onto the GUI
handles = copy_jtable(handles) ;

% --- Get the stats filename corresponding to selected data
function fname_stat = get_stats_fname(handles)

% Get selected data 
% Gets selected items in java table
all_idx = handles.table.jtable.getSelectedRows ; 
idx = [] ;
% Get proper indices in case table is sorted
for kk=1:length(all_idx); idx(kk) = str2num(handles.table.jtable.getValueAt(all_idx(kk),handles.table.jtable.getColumnCount - 1)) ; end

stats = [] ;

handles= update_data_table(handles);

% Get file name (data)
fname = cell2mat(handles.table.sFiles(idx(1))) ;

% Get file name (stats)
[PATHSTR,NAME,EXT] = fileparts(fname);
NAME = strrep(NAME,'_data','_stats');
fname_stat = fullfile(PATHSTR,strcat(NAME,EXT)) ;

% --- Executes on selection change in data table.
function edit_stats(hObject, eventdata, handles)

% Get stats filename
fname_stat= get_stats_fname(handles) ;

if exist(fname_stat,'file')
    % Load data file
    load(fname_stat);

    % For all stat configuration computed
    for ss=1:length(stats)

       str{ss} = sprintf('Stats %d : p<%-5.3f\t\tnboot:%-10.4d\t\t\tBaseline :[%-5.3f,%0.3f]\t\t\tResult: %0.3f',...
           ss,stats(ss).pthresh,...
           stats(ss).nboot,stats(ss).baseline(1),stats(ss).baseline(2),...
           stats(ss).threshdur);

    end
else 
    % No statistics were computed
    str = '';
end
    set(handles.list_stats,'Value',1);
    set(handles.list_stats,'String',str);
    set(handles.list_stats,'Max',length(str)); % Make it so you can select more than 1.
% end

% --- Executes on selection change in list_patient.
function list_patient_Callback(hObject, eventdata, handles)

if isfield(handles,'table') 
    
    if ~isempty(handles.table.mia_table) 

    persistent chk

    % Get all patient selected 
    idx_selected = get(hObject,'Value'); 
    list_patients = get(hObject,'String');

    selected_patients = list_patients(idx_selected);

    % No patient in the database
    if isempty(selected_patients) ; return ; end 
    
    % Highligths the data corresponding to patient 
    handles.table.jtable.clearSelection ; 

    % Get indices of row to select (if filters are used it is not
    % straitforward)
    for pp=1:length(selected_patients)

        % Get diplayed table 
        for jj=0:handles.table.jtable.getRowCount-1

            if strcmp(selected_patients(pp),handles.table.jtable.getValueAt(jj,0))
                handles.table.jtable.addRowSelectionInterval(jj,jj)
            end
        end
    end

    if isempty(chk)
        chk = 1;
        pause(0.4); %Add a delay to distinguish single click from a double click
        if chk == 1
              chk = []; % Simple click
        end
    else
        chk = [];
        % Create progress bar for patients processing
        hwait = waitbar(50,'','Units','Normalized','Name','Loading data....');
        % Make the waitbar stay on top
        set(hwait,'WindowStyle','modal')

        % Get patient selected 
        pt =get(handles.list_patient, 'Value') ;
        list_pt = get(handles.list_patient,'String') ;

        % Create file name 
        tmp = dir(fullfile(handles.extOPTIONS.outdir,list_pt{pt(1)},'*signal_LFP.mat')) ;
        fname = fullfile(handles.extOPTIONS.outdir,list_pt{pt(1)},tmp.name) ;

        % Call sanity GUI 
        mia_sanity_check_gui({fname});
     
        % Close progress bar
        delete(hwait) ;

    end

end

end

% --- Executes during object creation, after setting all properties.
function list_patient_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uipanel_workingdir_CreateFcn(hObject, eventdata, handles)

% --- Executes on button press in browse_database.
function browse_database_Callback(hObject, eventdata, handles)

% Open a directory browser
% directoryname = uigetdir(fullfile(fileparts(which(mfilename))), 'Pick a Directory');
dirname=get(handles.outdir,'String');

directoryname = uigetdir(dirname, 'Pick a Directory');

%use of cancel button => back to the next state
if directoryname==0
    set(handles.outdir,'String',dirname);
    return;
end

% Set the new working directory
handles = set_workingdirectory(handles,directoryname) ; 

guidata(hObject,handles);

% --- Set a new working directory
function handles = set_workingdirectory(handles,directoryname)

% Update MAINDIR text panel
set(handles.outdir,'String',directoryname);
handles.extOPTIONS.outdir = directoryname ;

handles = initialize_gui(handles.figure1,handles,false);

% Get User home directory
dirname = mia_getuserdir ;

handles.history.dirname = handles.extOPTIONS.outdir  ; 

history = handles.history ; 
save(fullfile(dirname,'.mia_history.mat'),'history');

% % --- Set a new loc table filename
% function handles = set_loctablefilename(handles,directoryname)
% 
% % Update MAINDIR text panel
% set(handles.outdir,'String',directoryname);
% handles.extOPTIONS.outdir = directoryname ;
% 
% handles = initialize_gui(handles.figure1,handles,false);
% 
% % Get User home directory
% dirname = getuserdir ;
% 
% handles.history.dirname = handles.extOPTIONS.outdir  ; 
% 
% history = handles.history ; 
% save(fullfile(dirname,'.mia_history.mat'),'history');


% Update the data list 
function handles = update_data_table(handles)

% Gets selected items selected in table
all_idx = handles.table.jtable.getSelectedRows ; 

% Get table that contains all files (names) found in the working directory
[handles.table.mia_table,handles.table.sFiles] = mia_create_table_workdir(handles.extOPTIONS.outdir, handles.current_loctable) ;

jtable = com.jidesoft.grid.SortableTable(handles.table.mia_table,{'Patient','Method','Montage','Freq. band','Fs','Remove Avg','Nb stats','Localized Contacts','ID'});

% These lines avoid Java Null Exception 
jtable.getColumnModel().getColumn(jtable.getColumnCount - 1).setMinWidth(0);
jtable.getColumnModel().getColumn(jtable.getColumnCount - 1).setMaxWidth(0);
jtable.getColumnModel().getColumn(jtable.getColumnCount - 1).setWidth(0);

handles.table.jtable.setModel(jtable.getModel()) ;

% Trick to hide to indexing column
handles.table.jtable.getColumnModel().getColumn(jtable.getColumnCount - 1).setMinWidth(0);
handles.table.jtable.getColumnModel().getColumn(jtable.getColumnCount - 1).setMaxWidth(0);
handles.table.jtable.getColumnModel().getColumn(jtable.getColumnCount - 1).setWidth(0);

% Retrieve selected items (Highligths)
if ~isempty(all_idx) 
    
     % Highligths the data corresponding to patient 
    handles.table.jtable.clearSelection ; 

    % Hightligth new entries in jtable
    % Get indices of row to select (if filters are used it is not
    % straitforward)
    for iFile=1:length(all_idx)

        % Get diplayed table 
        for iRow=0:handles.table.jtable.getRowCount-1
            if all_idx(iFile)==iRow
                handles.table.jtable.addRowSelectionInterval(iRow,iRow)
            end
        end
    end
end


% --- Executes on button press in display_sanity.
function display_sanity_Callback(hObject, eventdata, handles)

pt =get(handles.list_patient, 'Value') ;

% No file was selected 
if isempty(pt)
    errordlg('You must select a patient','Error');
else

    % Create progress bar for patients processing
    hwait = waitbar(50,'','Units','Normalized','Name','Loading data....');
    
    % Make the waitbar stay on top
    set(hwait,'WindowStyle','modal')

    % Get list of patients
    list_pt = get(handles.list_patient,'String') ;
    
    % Display all selected patients
    for iPt=1:length(pt)
        tmp = dir(fullfile(handles.extOPTIONS.outdir,list_pt{pt(iPt)},'*signal_LFP.mat')) ;
        fname = fullfile(handles.extOPTIONS.outdir,list_pt{pt(iPt)},tmp.name) ; 
        
        % Call main visualization GUI 
        % Conversion {} required (prevents warning : "The input to STR2FUNC...")
        mia_sanity_check_gui({fname});
        
    end
    
    % Close progress bar
    delete(hwait) ;

end  

% --- Executes on button press in display_stats.
function display_stats_Callback(hObject, eventdata, handles)

% ASD : is this line usefull?
[handles.table.mia_table,handles.table.sFiles] = mia_create_table_workdir(handles.extOPTIONS.outdir, handles.current_loctable) ;

% Gets selected items in java table
all_idx = handles.table.jtable.getSelectedRows ; 
idx = [] ;
% Get proper indices in case table is sorted
for kk=1:length(all_idx); idx(kk) = str2num(handles.table.jtable.getValueAt(all_idx(kk),handles.table.jtable.getColumnCount - 1)) ; end

% No file was selected 
if isempty(idx)
    errordlg('You must select a file','Error');
else

    % Loop throught all file selected in the table
    for ii=1:length(idx)
        % Create progress bar for patients processing
        hwait_pt = waitbar(0,'','Units','Normalized','Name','Loading data...');
        % Make the waitbar stay on top
        set(hwait_pt,'WindowStyle','modal')

        % Get selected patient file name 
        fname = handles.table.sFiles(idx(ii)) ;
        
        % Call main visualization GUI 
        % Conversion {} required (prevents warning : "The input to STR2FUNC...")
        mia_display_images_stats_gui(fname,handles.table.mia_table(idx(ii),:));

        % Load is done : close progress bar
        delete(hwait_pt) ;

    end
    
end

% --- Executes on button press in import_data.
function import_data_Callback(hObject, eventdata, handles)

% Error if no working dir
if ~isfield(handles.extOPTIONS,'outdir')
    warndlg('You must pick up a working directory first') ;
    return ;
end

% Import new patients data
handles = mia_import_data(handles) ;

% Update Patient table
handles = update_patientlist(handles) ;

guidata(hObject,handles);

% --- Update patient list
function handles = update_patientlist(handles) 

% Reads all folders that are in MAINDIR
d = dir(get(handles.outdir,'String'));
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

set(handles.list_patient,'Value',1);
set(handles.list_patient,'String',subjects);
set(handles.list_patient,'Max',length(subjects)); % Make it so you can select more than 1.
            
% --- Update study list
function handles = update_studies_list(handles) 

% Reads all folders that are in MAINDIR
group_foldnam = fullfile(fileparts(get(handles.outdir,'String')),'GA_Results') ; 

d = dir(group_foldnam);
isub = [d(:).isdir]; % returns logical vector if is folder
studies= {d(isub).name}';
studies(ismember(studies,{'.','..'})) = []; % Removes . and ..

set(handles.list_studies,'Value',1);
set(handles.list_studies,'String',studies);
% set(handles.list_study,'Max',length(studies)); % Make it so you can select more than 1.

% --- Update study list
function handles = update_atlas_list(handles) 

% Reads all folders that are in MAINDIR
group_foldnam = fullfile(fileparts(get(handles.outdir,'String')),'GA_Results') ; 

d = dir(group_foldnam);
isub = [d(:).isdir]; % returns logical vector if is folder
atlas= {d(~isub).name}';
atlas(ismember(atlas,{'.','..','.DS_Store'})) = []; % Removes systems files
atlas = strrep(atlas,'m_table_','') ; 
atlas = strrep(atlas,'.mat','') ; 

% If any localization table 
if ~isempty(atlas)
    set(handles.list_atlas,'Value',1);
    set(handles.list_atlas,'String',atlas);
    handles.current_loctable = load(fullfile(group_foldnam,strcat('m_table_',atlas{1})));
end
% set(handles.list_study,'Max',length(studies)); % Make it so you can select more than 1.

% --- Executes on button press in pushbutton_statistics.
function pushbutton_statistics_Callback(hObject, eventdata, handles)

handles = compute_statistics(handles) ;

% --- Executes on menu item "Statistics..." 
function menu_compute_stats_Callback(hObject, eventdata, handles)

handles = compute_statistics(handles) ;

% --- Compute statistics 
function [handles] = compute_statistics(handles) 

% Error if working directory was not set
if ~isfield(handles.extOPTIONS,'outdir')
    warndlg('You must pick up a working directory first') ;
    return ;
end

% Gets selected items in java table
idx = [] ;
if isfield(handles,'jtable')
    all_idx = handles.table.jtable.getSelectedRows ; 

    % Get proper indices in case table is sorted
    for kk=1:length(all_idx); idx(kk) = str2num(handles.table.jtable.getValueAt(all_idx(kk),handles.table.jtable.getColumnCount - 1)) ; end
else
    idx = 1;
end

if isempty(idx)
    errordlg('You must select a file','Error');
    return ;
end

% Call statistic analysis GUI
handles_stats = mia_statistics_gui(idx,handles.extOPTIONS.outdir,  handles.current_loctable);

% Update the list 
handles = update_patientlist(handles) ; %% ASD : WONT WORK??
handles= update_data_table(handles);

%--- Removes localisation table for a specific patient
function [handles] = resetPatientLoc(handles)

    % Get list of ALL patients
    list_patients = get(handles.list_patient,'String');

    % Delete all m_table.mat from patients directories
    for pp=1:length(list_patients)
        delete(fullfile(get(handles.outdir,'String'),list_patients{pp},'m_table.mat'));
    end            
   
    % Update data table
    handles= update_data_table(handles);

    % Update Patient table
    handles = update_patientlist(handles) ;


% % --- Executes on button press in pushbutton_statistics.
function pushbutton_loadloctable_Callback(hObject, eventdata, handles)


inputFormat{2,1} = {'xlsx'};
inputFormat{2,2} = 'Excel table (.xlsx)';

inputFormat{1,1} = {'tsv'};
inputFormat{1,2} = 'Brainstorm file (.tsv)';

% Open a Java file dialog box to browse a table 
[RawFile, FileFormat] = mia_dialog_getfile('MIA : Pick a labeling table...', ...  % Window title
                                            handles.extOPTIONS.outdir, ...          % Working directory
                                            inputFormat);    % List of available file formats

% Labeling table is in Excel format
if strcmp(FileFormat,'Excel table (.xlsx)')
    [struct_table, status, message] = mia_read_loc_table(RawFile{1}) ;

% Labeling table is in Braintorm .tsv format
elseif strcmp(FileFormat,'Brainstorm file (.tsv)')
    OPTIONS.patients = get(handles.list_patient,'String');
    if isempty(OPTIONS.patients)  
        errordlg('You must import patient data into MIA first','Error');
        return ;
    end
    [struct_table, status, message] = mia_read_loc_tsv_table(RawFile{1}, OPTIONS) ;

else ; return ; 

end
   
% Return error if doublons exist
if status==0
    message = sprintf('%s  :\n%s\n\n%s','Some contacts were found twice in the table : ',message,'Table was NOT loaded : Please remove duplicates and try again');
    errordlg(message) ;
    return
% Return error if no header in file
elseif status == -1 
    errordlg(message);   
    return
% User cancelled operation
elseif status == 2 
    return
else
    grpOPTIONS.maindir =  handles.extOPTIONS.outdir ; 
    % Map the contacts from loc table with the ones in the data 
    [s.m_table_all, status, message] = mia_get_dataloc_table(struct_table,grpOPTIONS);
end

if status==0
    warndlg(message) ;
   %return
end    
 % Prompt user for a study name
[~,NAME,~] = fileparts(RawFile{1}) ; 
def = {NAME,'hsv'};
loctable_name = mia_newid({'Enter an Atlas name:'},'ATLAS NAME',1,def);

% Use of cancel button : nothing happens
if isempty(loctable_name); return; end

resetPatientLoc(handles) ; 

% GA_Results directory 
group_dir=fullfile(fileparts(grpOPTIONS.maindir),'GA_Results') ; 

% Create GA_Results if does not exist
if ~exist(group_dir,'dir') ; mkdir(group_dir); end

% Atlas filename 
fname = cell2mat(fullfile(group_dir,strcat('m_table_',loctable_name,'.mat')));
s.table_fname = RawFile{1} ;

% Save m_table fopr this atlas
save(fname,'-struct','s');

handles = update_atlas_list(handles) ;


% --- Executes on button press in pushbutton_newStudy.
function pushbutton_newStudy_Callback(hObject, eventdata, handles)

if ~isfield(handles.extOPTIONS,'outdir')
    warndlg('You must pick up a working directory first') ;
    return ;
end

idx = [] ;

% ASD : DOES NOT WORK??

% Gets selected items in java table
if isfield(handles,'jtable')
    all_idx = handles.table.jtable.getSelectedRows ; 

    % Get proper indices in case table is sorted
    for kk=1:length(all_idx); idx(kk) = str2num(handles.table.jtable.getValueAt(all_idx(kk),handles.table.jtable.getColumnCount - 1)) ; end

else
    idx = 1;
end

if isempty(idx)
    errordlg('You must select a value','Error');
    return ;
end

% Call create_a_study GUI (guide)
mia_create_a_study_gui(idx,handles.extOPTIONS.outdir, handles.current_loctable);

% Update jtable
handles= update_data_table(handles);
handles = update_studies_list(handles) ;

guidata(hObject, handles);

% --- Executes on button press in GO.
function GO_Callback(hObject, eventdata, handles)

[PATH,~,~]=fileparts(handles.extOPTIONS.outdir);
maindir=char(fullfile(PATH,'GA_Results'));

if ~exist(maindir,'dir')
    warndlg('You must create a study first') ;
    return ;
else
    d = dir(maindir);
    study = {d.name}';
    study(ismember(study,{'.','..','.DS_Store'})) = []; %remove . et ..
    if isempty(study)
        warndlg('You must create a study first') ;
        return ;
    end
end

% Read selected labelling Atlas
idx_selected = get(handles.list_atlas,'Value'); 

% Get selected labelling atlas 
list_altas = get(handles.list_atlas,'String');

% There is no m_table (no Atlas) 
if isempty(list_altas)
     warndlg('You must create an atlas first') ;
        return ;
end


selected_atlas = list_altas(idx_selected);

load(cell2mat(fullfile(maindir, strcat('m_table_',selected_atlas))));

% ASD : add message expliciting the atlas that was taken for group GUI

% Start group GUI
mia_ganalysis_gui(m_table_all,maindir,selected_atlas, handles.extOPTIONS);


% --- Executes when figure1 is resized.
function [handles] = figure1_SizeChangedFcn(hObject, eventdata, handles)

% Update jtable size only if it exists
if isfield(handles,'table');
    if ~isempty(handles.table.mia_table);
        [handles] = copy_jtable(handles);
       
    end
end

% Copy the java table on the main panel
function [handles] = copy_jtable(handles)

 % Temporarly sets panel units to Pixels
set(handles.uipanel_listdata,'Units','Pixels');

% Get the position in Pixels 
y = get(handles.uipanel_listdata,'Position');

% Fits the jtable on the panel
[handles.hjtable,handles.hjcontainer] = javacomponent(handles.table.jScrollPane,[0,0,y(3),y(4)],handles.uipanel_listdata);

% Sets back the panel in normalized units
set(handles.uipanel_listdata,'Units','normalized');

% --- Executes on button press in pushbutton_plusImport (import data)
function pushbutton_plusImport_Callback(hObject, eventdata, handles)

% Error if no working dir
if ~isfield(handles.extOPTIONS,'outdir')
    warndlg('You must pick up a working directory first') ;
    return ;
end

handles = mia_import_data(handles) ;

% Update Patient table
handles = update_patientlist(handles) ;


% --- Executes on button press in pushbutton_minusRemovePt.
function pushbutton_minusRemovePt_Callback(hObject, eventdata, handles)

% Get all patient selected 
idx_selected = get(handles.list_patient,'Value'); 
list_patients = get(handles.list_patient,'String');

selected_patients = list_patients(idx_selected);

str = sprintf('Are you sure you want to remove %s patient(s)?',...
      num2str(length(idx_selected)));

    ButtonName = questdlg(str, ...
        'Remove patient','Yes', 'No','No');
    switch ButtonName,
        case 'Yes',
           for pp=1:length(selected_patients)
             rmdir(char(fullfile(handles.extOPTIONS.outdir,selected_patients(pp))),'s');
            end            
        case 'No'
            return;
           
    end

% Update jtable
handles= update_data_table(handles);

% Update Patient table
handles = update_patientlist(handles) ;

guidata(hObject,handles);

% --- Execute on pushbutton extract frequency.
function menu_freq_extract_Callback(hObject, eventdata, handles)

handles = compute_frequencies(handles) ;

% --- Executes on button press in pushbutton_plusExtractFreq.
function pushbutton_plusExtractFreq_Callback(hObject, eventdata, handles)

handles = compute_frequencies(handles) ;

% --- Executes on button press in pushbutton_ExtractFreq.
function pushbutton_ExtractFreq_Callback(hObject, eventdata, handles)

handles = compute_frequencies(handles) ;

% Update the table
handles = update_data_table(handles);

% --- Extract frequencies and highlight results in jtable
function [handles] = compute_frequencies(handles)
% Error if no working dir
if ~isfield(handles.extOPTIONS,'outdir')
    warndlg('You must pick up a working directory first') ;
    return ;
end

% Call frequency analysis GUI
newEntries= mia_extract_frequency(handles.list_patient, handles.extOPTIONS.outdir);

% Update the table
handles = update_data_table(handles);

% If any new element of the table
if ~isempty(newEntries) 
    
     % Highligths the data corresponding to patient 
    handles.table.jtable.clearSelection ; 

    % Hightligth new entries in jtable
    % Get indices of row to select (if filters are used it is not
    % straitforward)
    for iFile=1:length(newEntries)

        % Get diplayed table 
        for iRow=0:handles.table.jtable.getRowCount-1
            iFileInTab = str2num(handles.table.jtable.getValueAt(iRow,handles.table.jtable.getColumnCount - 1)) ; 
            if strcmp(newEntries(iFile),handles.table.sFiles(iFileInTab))
                handles.table.jtable.addRowSelectionInterval(iRow,iRow)
            end
        end
    end
end

% --- Executes on button press in pushbutton_minusRmData.
function pushbutton_minusRmData_Callback(hObject, eventdata, handles)

[handles.table.mia_table,handles.table.sFiles] = mia_create_table_workdir(handles.extOPTIONS.outdir, handles.current_loctable) ;

all_idx = handles.table.jtable.getSelectedRows ;

idx = [] ;
% Get proper indices in case table is sorted
for kk=1:length(all_idx); idx(kk) = str2num(handles.table.jtable.getValueAt(all_idx(kk),handles.table.jtable.getColumnCount - 1)) ; end

% No file was selected 
if isempty(idx)
    errordlg('You must select a file','Error');
else

    str = sprintf('Are you sure you want to remove %s data file(s)?',...
          num2str(length(idx)));

    ButtonName = questdlg(str, ...
            'Remove patient','Yes', 'No','No');
    switch ButtonName,
        case 'Yes',
           for pp=1:length(idx)
               fprintf(sprintf('Remove %s\n',char(handles.table.sFiles(idx(pp)))));
               delete(char(handles.table.sFiles(idx(pp))));
               
           end            
        case 'No'
            return;
           
    end

    % Update jtable
    handles= update_data_table(handles);
    
    % Update Patient table
    handles = update_patientlist(handles) ;

    guidata(hObject,handles);

end

% --------------------------------------------------------------------
function about_Callback(hObject, eventdata, handles)
mia_about(handles) ;

% --------------------------------------------------------------------
function menu_import_patient_Callback(hObject, eventdata, handles)
% hObject    handle to menu_import_patient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Error if no working dir
if ~isfield(handles.extOPTIONS,'outdir')
    warndlg('You must pick up a working directory first') ;
    return ;
end

% Import new patients data
handles = mia_import_data(handles) ;

% Update Patient table
handles = update_patientlist(handles) ;

guidata(hObject,handles);

% --------------------------------------------------------------------
function menu_preferences_Callback(hObject, eventdata, handles)
% hObject    handle to menu_preferences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  warndlg('This menu is not implemented yet.') ;

% --- Executes on button press in pushbutton_minusRmLocTable.
function pushbutton_minusRmLocTable_Callback(hObject, eventdata, handles)

% Get patient selected 
idx_selected =get(handles.list_patient, 'Value') ;
list_pt = get(handles.list_patient,'String') ;

selected_patients = list_pt(idx_selected);

if isempty(selected_patients) 
    errordlg('You must select a patient','Error');
 
else

    spt = sprintf('%s\n',selected_patients{:}); 
    str = sprintf('%s%s\n','Are you sure you want to remove localizations file(s) for ',spt);

    ButtonName = questdlg(str, ...
            'Remove localizations','Yes', 'No','No');
    switch ButtonName,
        case 'Yes',
           for pp=1:length(selected_patients)
             delete(char(fullfile(get(handles.outdir,'String'),selected_patients(pp),'m_table.mat')));
           end            
        case 'No'
            return;

    end

    % Update jtable
    handles= update_data_table(handles);
 
    % Update Patient table
    handles = update_patientlist(handles) ;

    guidata(hObject,handles);
end 


% --------------------------------------------------------------------
function menu_quit_Callback(hObject, eventdata, handles)
close(handles.figure1);


% --------------------------------------------------------------------
function menu_help_Callback(hObject, eventdata, handles)


% --- Executes on selection change in list_studies.
function list_studies_Callback(hObject, eventdata, handles)
% hObject    handle to list_studies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_studies contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_studies


% --- Executes during object creation, after setting all properties.
function list_studies_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_studies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in list_atlas.
function list_atlas_Callback(hObject, eventdata, handles)

persistent chk

if isfield(handles,'table') 
    
    if ~isempty(handles.table.mia_table) && ~isempty(handles.current_loctable)

        % Get table selected 
        idx_selected = get(hObject,'Value'); 
        list_table = get(hObject,'String');

        selected_table = list_table(idx_selected);

        % Reads all folders that are in MAINDIR
        group_foldnam = fullfile(fileparts(get(handles.outdir,'String')),'GA_Results') ; 

        d = dir(group_foldnam);
        isub = [d(:).isdir]; % returns logical vector if is folder
        loctables= {d(~isub).name}';
        loctables(ismember(loctables,{'.','..','.DS_Store'})) = []; % Removes systems files
        loctables = strrep(loctables,'m_table_','') ; 
        loctables = strrep(loctables,'.mat','') ; 

        handles.current_loctable = load(fullfile(group_foldnam,strcat('m_table_',loctables{idx_selected})));

        handles= update_data_table(handles);
                        
         if isempty(chk)
            chk = 1;
            pause(0.2); %Add a delay to distinguish single click from a double click
            if chk == 1
                  chk = []; % Simple click

            end
        else
            chk = [];
             % Prompt user for a study name
            prompt = {'Enter an ROIs table name:'};
            dlg_title = 'ROI TABLE NAME';
            num_lines = 1;
            def = {strrep(selected_table{1},'.xlsx',''),'hsv'};
            loctable_name = newid(prompt,dlg_title,num_lines,def);

            % Atlas filename 
            fname_old = cell2mat(fullfile(group_foldnam,strcat('m_table_',selected_table,'.mat')));
            fname_new = cell2mat(fullfile(group_foldnam,strcat('m_table_',loctable_name,'.mat')));

            % Use of cancel button OR use same name as previous one
            if isempty(loctable_name)||strcmp(fname_old,fname_new)==1
                return;
            end

            % Save m_table fopr this atlas
            movefile(fname_old,fname_new);

            handles = update_atlas_list(handles) ;

         end
         % Update text field with localization table path
         handles = update_loctablepath_text(handles);

    end
    
end
   
    
function handles = update_loctablepath_text(handles) 
% Displays path of localization table 
% Not saved in previous MIA version, if unknown write "unknown"
if isfield(handles.current_loctable,'table_fname')
    set(handles.text_labellingFile,'String',handles.current_loctable.table_fname);
else
    set(handles.text_labellingFile,'String','Unknown');
end
        
% --- Executes during object creation, after setting all properties.
function list_atlas_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_remove_atlas.
function pushbutton_remove_atlas_Callback(hObject, eventdata, handles)

% Get patient selected 
idx_selected =get(handles.list_atlas, 'Value') ;
list_atlas = get(handles.list_atlas,'String') ;
selected_atlas = list_atlas(idx_selected);

[PATH,~,~]=fileparts(handles.extOPTIONS.outdir);
group_dir=char(fullfile(PATH,'GA_Results'));

if isempty(selected_atlas) 
    errordlg('You must select an atlas','Error');
else
    
    str = sprintf('Are you sure you want to remove %s ?',cell2mat(selected_atlas));
    
    button = questdlg(str, 'Remove Atlas','Yes', 'No','No');
    switch button,
        case 'Yes',
             delete(char(fullfile(group_dir,strcat('m_table_',selected_atlas,'.mat'))));
          
        case 'No'
            return;

    end

    % Update Atlas table
    handles = update_atlas_list(handles) ;

    guidata(hObject,handles);

end

% --- Executes on button press in pushbutton_remove_study.
function pushbutton_remove_study_Callback(hObject, eventdata, handles)

% Get all patient selected 
idx_selected = get(handles.list_studies,'Value'); 
list_studies = get(handles.list_studies,'String');
selected_studies = list_studies(idx_selected);

% Get studies path
[PATH,~,~]=fileparts(handles.extOPTIONS.outdir);
group_dir=char(fullfile(PATH,'GA_Results'));

str = sprintf('Are you sure you want to remove %s ?',cell2mat(selected_studies)); 
button = questdlg(str, 'Remove Study','Yes', 'No','No');

switch button,
    case 'Yes',
         rmdir(char(fullfile(group_dir,selected_studies)),'s');
    case 'No'
        return;

end

% Update jtable
handles= update_studies_list(handles);

guidata(hObject,handles);


% --- Executes on selection change in list_stats.
function list_stats_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function list_stats_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_remove_stats.
function pushbutton_remove_stats_Callback(hObject, eventdata, handles)

% Get all patient selected 
idx_selected = get(handles.list_stats,'Value'); 

str = sprintf('Are you sure you want to remove this stat ?'); 
button = questdlg(str, 'Remove Stats','Yes', 'No','No');


switch button,
    case 'Yes',
         % Get stats filename
        fname_stat= get_stats_fname(handles) ;

        if exist(fname_stat,'file')

            % Load data file
            load(fname_stat);
            stats(idx_selected) = [] ; 
            save(fname_stat,'stats');
            edit_stats(hObject, eventdata, handles);
            update_data_table(handles);
        end 

    case 'No'
        return;

end


% --- Executes during object creation, after setting all properties.
function pushbutton_loadloctable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadloctable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function MIA_website_Callback(hObject, eventdata, handles)
% hObject    handle to MIA_website (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 web('http://www.neurotrack.fr/mia/', '-browser')
 


% --------------------------------------------------------------------
function troubleshoot_Callback(hObject, eventdata, handles)
% hObject    handle to troubleshoot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web("mailto:anne-sophie.dubarry@univ-amu.fr") 
