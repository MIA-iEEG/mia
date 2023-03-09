function varargout = mia_ganalysis_gui(varargin)
% MIA_GANALYSIS_GUI MATLAB code for mia_ganalysis_gui.fig
%      MIA_GANALYSIS_GUI, by itself, creates a new MIA_GANALYSIS_GUI or raises the existing
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
% Copyright (C) 2016-2022 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)

gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @mia_ganalysis_gui_OpeningFcn, ...
    'gui_OutputFcn',  @mia_ganalysis_gui_OutputFcn, ...
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


% --- Executes just before mia_ganalysis_gui is made visible.
function mia_ganalysis_gui_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for mia_ganalysis_gui
handles.output = hObject;

%handles.group_files=varargin{1};
handles.m_table_as = varargin{1};
handles.group_dir = varargin{2};
selected_atlas = varargin{3};
handles.db_dir = varargin{4};
handles.list_study_mia = varargin{5};

% Jtable position offset in Y (down)
handles.shift = 42;
handles.table_size = [15,10,735,600];

% Set distinct colors for pt
handles.dOPTIONS.clr = jet(numel(unique(handles.m_table_as(:,1)))); 
       
% Set title with the name of the input atlas 
set(handles.text_title,'string',strcat('GROUP ANALYSIS : ',selected_atlas));

% Find all files containg group analysis (rois, etc.)
group_folders = fullfile(handles.group_dir,handles.list_study_mia.String,strcat(handles.list_study_mia.String,'*.mat')); 
for ii=1:length(group_folders)
    d = dir(group_folders{ii});
    handles.group_files{ii} = d.name;   
end

 % Here extract the name of the method out of handles.group_files(Study_nbpt_name)

% Update listbox
set(handles.list_study,'String',handles.list_study_mia.String);
set(handles.list_study,'Value',handles.list_study_mia.Value);

% Get current group file
current_grp_file = fullfile(handles.group_dir, handles.list_study_mia.String{handles.list_study_mia.Value}, handles.group_files{handles.list_study_mia.Value}) ; 

OPTIONS.allow_flipsign = get(handles.allow_flipsign,'value') ; 
OPTIONS.signifmode = get(handles.checkbox_significant_only,'value') ; 
OPTIONS.flip_thresh = str2double(get(handles.flip_thresh,'String')) ; 
   
 [handles.rois,...
     handles.table.mia_table,...
     handles.edges,handles.datafiles,...
     handles.mtg,...
     handles.freqb,...
     handles.wdir] = mia_create_table_of_rois(handles.m_table_as, current_grp_file, OPTIONS) ;
 
% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = mia_ganalysis_gui_OutputFcn(hObject, eventdata, handles)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in list_study.
function list_study_Callback(hObject, eventdata, handles)

%update jtable when study selected study get changed
if get(handles.togglebutton1,'Value')==1
    
    delete(handles.hjtable);
    guidata(hObject, handles);
    
    % Get current group file
    current_grp_file = fullfile(handles.group_dir, handles.list_study.String{handles.list_study.Value}, handles.group_files{handles.list_study.Value}) ; 
    
    OPTIONS.signifmode = get(handles.checkbox_significant_only,'value') ;
    OPTIONS.allow_flipsign = get(handles.allow_flipsign,'value') ; 
    OPTIONS.flip_thresh = str2double(get(handles.flip_thresh,'String')) ; 
   
    % Get table that contains all files (names) found in the working directory
    [handles.rois,...
     handles.table.mia_table,...
     handles.edges,handles.datafiles,...
     handles.mtg,...
     handles.freqb,...
     handles.wdir] = mia_create_table_of_rois(handles.m_table_as, current_grp_file,OPTIONS) ;
 
%     jtable = com.jidesoft.grid.SortableTable(handles.table.mia_table,{'Region','Onset','Patients Correlation','Channels Correlation','N patients','N contacts','ID'});
    jtable = com.jidesoft.grid.SortableTable(handles.table.mia_table,{'Region','Patients Correlation','Channels Correlation','N patients','N contacts','ID'});
    tableHeader = com.jidesoft.grid.AutoFilterTableHeader(jtable);
        
    % Trick to hide to indexing column
    jtable.getColumnModel().getColumn(jtable.getColumnCount - 1).setMinWidth(0);
    jtable.getColumnModel().getColumn(jtable.getColumnCount - 1).setMaxWidth(0);
    jtable.getColumnModel().getColumn(jtable.getColumnCount - 1).setWidth(0);

    tableHeader.setAutoFilterEnabled(true)
    tableHeader.setShowFilterName(true)
    tableHeader.setShowFilterIcon(true)
    jtable.setTableHeader(tableHeader)
    installer = com.jidesoft.grid.TableHeaderPopupMenuInstaller(jtable);
    pmCustomizer1=com.jidesoft.grid.AutoResizePopupMenuCustomizer;
    installer.addTableHeaderPopupMenuCustomizer(pmCustomizer1);
    pmCustomizer2=com.jidesoft.grid.TableColumnChooserPopupMenuCustomizer;
    installer.addTableHeaderPopupMenuCustomizer(pmCustomizer2);
    jScrollPane = javax.swing.JScrollPane(jtable);
    
    handles.table.jtable = jtable ;
    handles.table.jScrollPane = jScrollPane ;
    
    handles.jtable = jtable ;
    [handles.hjtable,handles.hjcontainer]=javacomponent(jScrollPane,handles.table_size,gcf);

else
    
    % Get current group file
    current_grp_file = fullfile(handles.group_dir, handles.list_study.String{handles.list_study.Value}, handles.group_files{handles.list_study.Value}) ; 

    OPTIONS.signifmode = get(handles.checkbox_significant_only,'value') ; 
    OPTIONS.allow_flipsign = get(handles.allow_flipsign,'value') ; 
    OPTIONS.flip_thresh = str2double(get(handles.flip_thresh,'String')) ; 
   
   [handles.rois,...
     handles.table.mia_table,...
     handles.edges,handles.datafiles,...
     handles.mtg,...
     handles.freqb,...
     handles.wdir] = mia_create_table_of_rois(handles.m_table_as, current_grp_file, OPTIONS) ;
end
  guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function list_study_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)


if (get(hObject,'Value') == get(hObject,'Max'))
    
    % Get current group file
    current_grp_file = fullfile(handles.group_dir, handles.list_study_mia.String{handles.list_study_mia.Value}, handles.group_files{handles.list_study_mia.Value}) ; 
    
    OPTIONS.signifmode = get(handles.checkbox_significant_only,'value') ; 
    OPTIONS.allow_flipsign = get(handles.allow_flipsign,'value') ; 
    OPTIONS.flip_thresh = str2double(get(handles.flip_thresh,'String')) ; 
    
    % Get table that contains all files (names) found in the working directory
      [handles.rois,...
     handles.table.mia_table,...
     handles.edges,handles.datafiles,...
     handles.mtg,...
     handles.freqb,...
     handles.wdir] = mia_create_table_of_rois(handles.m_table_as, current_grp_file,OPTIONS) ;
    
    % Enlarge main figure in height and shift down
    fig_position = get(handles.figure1,'Position') ;
    fig_position(4) = fig_position(4) + handles.shift ;
    fig_position(2) = fig_position(2) - handles.shift ;
    set(handles.figure1,'Position',fig_position) ;

    uipanel_displayrois_pos = get(handles.uipanel_displayrois,'Position');
    uipanel_displayrois_pos(2) = uipanel_displayrois_pos(2) + handles.shift;
    set(handles.uipanel_displayrois,'Position',uipanel_displayrois_pos);

    uipanel_study_pos = get(handles.uipanel_studies,'Position');
    uipanel_study_pos (2) = uipanel_study_pos (2) + handles.shift;
    set(handles.uipanel_studies,'Position',uipanel_study_pos );

    uipanel_text_pos = get(handles.text_title,'Position');
    uipanel_text_pos  (2) = uipanel_text_pos  (2) + handles.shift;
    set(handles.text_title,'Position',uipanel_text_pos  );

    uipanel_flip_pos = get(handles.uipanel_flip,'Position');
    uipanel_flip_pos(2) = uipanel_flip_pos(2) + handles.shift;
    set(handles.uipanel_flip,'Position',uipanel_flip_pos );

    uipanel_rasterplot_pos = get(handles.uipanel_rasterplot,'Position');
    uipanel_rasterplot_pos(2) = uipanel_rasterplot_pos(2) + handles.shift;
    set(handles.uipanel_rasterplot,'Position',uipanel_rasterplot_pos );

    togglebutton_pos = get(handles.togglebutton1,'Position');
    togglebutton_pos(2) = togglebutton_pos(2) + handles.shift;
    set(handles.togglebutton1,'Position',togglebutton_pos);

%     jtable = com.jidesoft.grid.SortableTable(handles.table.mia_table,{'Region','Onset','Patients Correlation','Channels Correlation','N patients','N contacts','ID'});
    jtable = com.jidesoft.grid.SortableTable(handles.table.mia_table,{'Region','Patients Correlation','Channels Correlation','N patients','N contacts','ID'});

    % Trick to hide to indexing column
    jtable.getColumnModel().getColumn(jtable.getColumnCount - 1).setMinWidth(0);
    jtable.getColumnModel().getColumn(jtable.getColumnCount - 1).setMaxWidth(0);
    jtable.getColumnModel().getColumn(jtable.getColumnCount - 1).setWidth(0);

    if ~isempty(handles.rois)
  
        tableHeader = com.jidesoft.grid.AutoFilterTableHeader(jtable);
        tableHeader.setAutoFilterEnabled(true)
        tableHeader.setShowFilterName(true)
        tableHeader.setShowFilterIcon(true)
        jtable.setTableHeader(tableHeader)
        installer = com.jidesoft.grid.TableHeaderPopupMenuInstaller(jtable);
        pmCustomizer1=com.jidesoft.grid.AutoResizePopupMenuCustomizer;
        installer.addTableHeaderPopupMenuCustomizer(pmCustomizer1);
        pmCustomizer2=com.jidesoft.grid.TableColumnChooserPopupMenuCustomizer;
        installer.addTableHeaderPopupMenuCustomizer(pmCustomizer2);
 
        
    else 
         fprintf('No ROI to show\n');
    end
    
    jScrollPane = javax.swing.JScrollPane(jtable);
    handles.table.jtable = jtable ;
    handles.table.jScrollPane = jScrollPane ;
        
    handles.jtable = jtable ;
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

    [handles.hjtable,handles.hjcontainer]=javacomponent(jScrollPane,handles.table_size,gcf);

    guidata(hObject, handles);

    
else
    
    % Delete Java table
    delete(handles.hjtable);
    guidata(hObject, handles);
    
    % Reduce main figure in height and shift up figure + all components
    fig_position = get(handles.figure1,'Position') ;
    fig_position(4) = fig_position(4) - handles.shift ;
    fig_position(2) = fig_position(2) + handles.shift ;
    set(handles.figure1,'Position',fig_position) ;
    
    uipanel_displayrois_pos = get(handles.uipanel_displayrois,'Position');
    uipanel_displayrois_pos(2) = uipanel_displayrois_pos(2) - handles.shift;
    set(handles.uipanel_displayrois,'Position',uipanel_displayrois_pos);
    
    uipanel_study_pos = get(handles.uipanel_studies,'Position');
    uipanel_study_pos (2) = uipanel_study_pos (2) - handles.shift;
    set(handles.uipanel_studies,'Position',uipanel_study_pos );
    
    uipanel_text_pos = get(handles.text_title,'Position');
    uipanel_text_pos  (2) = uipanel_text_pos  (2) -handles.shift;
    set(handles.text_title,'Position',uipanel_text_pos  );
    
    uipanel_rasterplot_pos = get(handles.uipanel_rasterplot,'Position');
    uipanel_rasterplot_pos(2) = uipanel_rasterplot_pos(2) - handles.shift;
    set(handles.uipanel_rasterplot,'Position',uipanel_rasterplot_pos );
    
    uipanel_flip_pos = get(handles.uipanel_flip,'Position');
    uipanel_flip_pos(2) = uipanel_flip_pos(2) - handles.shift;
    set(handles.uipanel_flip,'Position',uipanel_flip_pos );
 
    togglebutton_pos = get(handles.togglebutton1,'Position');
    togglebutton_pos(2) = togglebutton_pos(2) - handles.shift;
    set(handles.togglebutton1,'Position',togglebutton_pos);
    
    
end
  guidata(hObject, handles);


% --- Executes on button press in display_rois.
function display_rois_Callback(hObject, eventdata, handles)

if get(handles.togglebutton1,'Value')==0
    errordlg('Use show table button to select a region to explore','Error');
else
    % Gets selected items in java table
    all_idx = handles.table.jtable.getSelectedRows ; 
    idx = [] ;
    % Get proper indices in case table is sorted
    for kk=1:length(all_idx); idx(kk) = str2num(handles.table.jtable.getValueAt(all_idx(kk),handles.table.jtable.getColumnCount - 1)) ; end
    if isempty(idx)
        errordlg('You must select a region to explore','Error');
    % Several items selected in the table
    elseif length(idx)>0
  
        dOPTIONS.clr = jet(numel(unique(handles.m_table_as(:,1)))); % distinct colors for pt
        dOPTIONS.win_noedges = handles.edges;
        mia_display_roi(handles.rois(idx),dOPTIONS);
    end
end
  guidata(hObject, handles);


% --- Executes on button press in diplay_summary.
function diplay_summary_Callback(hObject, eventdata, handles)

% Pop a question window
ButtonName = questdlg('Display type : ', ...
    'Display ROI summary', ...
    'Colored time series','Stack plot','Flat display','Colored time series'); %,'Colored Time series with Ruban'

% No button was pressed (window closed)
if isempty(ButtonName); return ; end

if get(handles.togglebutton1,'Value')==0
    errordlg('Use show table button to select a region to explore','Error');
else
    % Gets selected items in java table
    all_idx = handles.table.jtable.getSelectedRows ; 
    idx = [] ;
    % Get proper indices in case table is sorted
    for kk=1:length(all_idx); idx(kk) = str2num(handles.table.jtable.getValueAt(all_idx(kk),handles.table.jtable.getColumnCount - 1)) ; end

      if isempty(idx)
        errordlg('You must select a region to explore','Error');
        % Several items selected in the table
    elseif length(idx)>0

    switch ButtonName,
       case 'Colored time series',
            pOPTIONS.thresh =1.96; % -1 for no color chronological organization
            pOPTIONS.threshdisp = 1.96; 
            pOPTIONS.nsub =1; % number of subplot
            pOPTIONS.title = '';
            [labels_o, colorm] = mia_display_summary_roi(handles.rois(idx),pOPTIONS) ;
        case 'Flat display',
            %pOPTIONS.thresh =1.96; % -1 for no color chronological organization
            pOPTIONS.thresh =3; % -1 for no color chronological organization
            pOPTIONS.threshdisp = 1.96; 
            pOPTIONS.nsub =1; % number of subplot
            pOPTIONS.title = '';
       
            [labels_o, colorm] = mia_display_summary_roi_flat(handles.rois(idx),pOPTIONS) ;
         
      case 'Stack plot',
           pOPTIONS.thresh =1.96; % -1 for no color chronological organization
            pOPTIONS.threshdisp = 1.96; 
            pOPTIONS.nsub =1; % number of subplot
           pOPTIONS.title = '';
       
           [labels_o, colorm] = mia_display_summary_roi_stack(handles.rois(idx),pOPTIONS) ;
          
%            while(1) 
                [x,y]=ginput(1); % get a click
                xlim = get(gca,'xlim');
                ylim = get(gca,'ylim');
                plot([x x], ylim); % vertical line
%            end
    end % end Switch    
      end % End if 
end
  guidata(hObject, handles);


% --- Executes on button press in one_patient_pushbutton.
function one_patient_pushbutton_Callback(hObject, eventdata, handles)

if get(handles.togglebutton1,'Value')==0
    errordlg('Use show table button to select a region to explore','Error');
else
    % Gets selected items in java table
    all_idx = handles.table.jtable.getSelectedRows ; 
    idx = [] ;
    % Get proper indices in case table is sorted
    for kk=1:length(all_idx); idx(kk) = str2num(handles.table.jtable.getValueAt(all_idx(kk),handles.table.jtable.getColumnCount - 1)) ; end

    if isempty(idx)
        errordlg('You must select a region to explore','Error');
        % Several items selected in the table
    elseif length(idx)>0
        
        froi=handles.rois(idx);
        f = [froi{:}];
        pt_votes = vertcat(f.namePt) ;
        
        %Gets all unique occurence of patients
        un = unique(pt_votes);
        
        [s,~]=listdlg('PromptString','Select a patient :', 'SelectionMode','multiple','ListString',un,'Name','One Patient Data');
        
        if ~isempty(s)
            un=un(s);
            for pp=1:length(un)
                
                pOPTIONS.ptKey = un{pp};
                pOPTIONS.title = pOPTIONS.ptKey;
                pOPTIONS.thresh = 3; % -1 for no color chronological organization
                pOPTIONS.nsub =1; % number of subplot
                
                % Display data for one specific patient
                [labels_o_pt, colorm_pt] =mia_display_summary_roi_pt(froi,pOPTIONS) ;
                
            end
        end
    end
end
guidata(hObject, handles);


% --- Executes on button press in TRfiles.
function TRfiles_Callback(hObject, eventdata, handles) %#ok<*DEFNU>

hwarn = warndlg('This will erase previous RT (Press cancel on next window to cancel)');
waitfor(hwarn);

% Open a directory browser
[filename, pathname] = uigetfile('*.xlsx', 'Pick a EXCEL file',handles.db_dir);

% CANCEL
if filename==0 ; return; end 

% % FOR DEBUG ONLY
% filename = '20200124_ BehavData_Dorsal_reformattedMIA.xlsx' ; 
% pathname =  '/Users/anne-sophiedubarry/Documents/2_DATA/Dorsal/' ; 
    
% Read the RT table
[rt_byPatients] = mia_read_rt_table(fullfile(pathname,filename));

% Get montage (from file name) 
rOPTIONS.freq= handles.freqb ;
rOPTIONS.mtg = handles.mtg ;
rOPTIONS.maindir = handles.db_dir;%handles.wdir; % ASD 2018/6/25 : TODO : Remove handles.wdir elsewhere (no use) 
rOPTIONS.freq = strcat(num2str(handles.freqb(1)),'_',num2str(handles.freqb(end)));
rOPTIONS.win_noedges = handles.edges;
rOPTIONS.clr = handles.dOPTIONS.clr ;

% Inject the RTs (organized by patients at this point) into ROIs
[handles.rois] = mia_add_rts(handles.rois,rt_byPatients, rOPTIONS) ; 

% TODO : ASD Here save new ROI strucutre with Rts
% for next exceution (if not we will have to reload RT
% at every execution) 

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in load_rasters.
function load_rasters_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to load_rasters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get montage (from file name) 
rOPTIONS.mtg = handles.mtg ;
%rOPTIONS.freq= handles.freqb ;
rOPTIONS.maindir = handles.db_dir;%handles.wdir; % ASD 2018/6/25 : TODO : Remove handles.wdir elsewhere (no use) 
rOPTIONS.freq = strcat(num2str(handles.freqb(1)),'_',num2str(handles.freqb(end)));
rOPTIONS.win_noedges = handles.edges;
rOPTIONS.clr = handles.dOPTIONS.clr ;

% if isempty(handles.method) ; rOPTIONS.method = 'hilbert' ; else rOPTIONS.method = handles.method ; end

% Get signals (single trials) 
[handles.rois] = mia_get_rasters(handles.rois,rOPTIONS) ;

% TODO : ASD Here save new ROI strucutre with Rts
% for next exceution (if not we will have to reload signals
% at every execution) 

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in display_rasters.
function display_rasters_Callback(hObject, eventdata, handles)

if get(handles.togglebutton1,'Value')==0
    errordlg('Use show table button to select a region to explore','Error');
else
    % Gets selected items in java table
    all_idx = handles.table.jtable.getSelectedRows ; 
    idx = [] ;
    % Get proper indices in case table is sorted
    for kk=1:length(all_idx); idx(kk) = str2num(handles.table.jtable.getValueAt(all_idx(kk),handles.table.jtable.getColumnCount - 1)) ; end
    if isempty(idx)
        errordlg('You must select a region to explore','Error');
    % One or Several items selected in the table
    elseif length(idx)>0
        % Get montage (from file name) 
        rOPTIONS.mtg = handles.mtg ;
        rOPTIONS.freq= handles.freqb ;
        rOPTIONS.maindir = handles.db_dir;%handles.wdir; % ASD 2018/6/25 : TODO : Remove handles.wdir elsewhere (no use) 
        rOPTIONS.freq = strcat(num2str(handles.freqb(1)),'_',num2str(handles.freqb(end)));
        rOPTIONS.win_noedges = handles.edges;
        rOPTIONS.clr = handles.dOPTIONS.clr ;

        mia_display_roi_and_raster(handles.rois(idx),rOPTIONS) ;
    end
end

 
% --- Executes during object creation, after setting all properties.
function uipanel_displayrois_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function uipanel_rasterplot_CreateFcn(hObject, eventdata, handles)

% --- Executes on button press in allow_flipsign.
function allow_flipsign_Callback(hObject, eventdata, handles)

% Hide (toggle callback) table if check button was checked or unchecked
if  get(handles.togglebutton1,'Value')== 1
    set(handles.togglebutton1,'Value',0); 
    togglebutton1_Callback(handles.togglebutton1, eventdata, handles);
end

function flip_thresh_Callback(hObject, eventdata, handles)

flip_thresh = str2double(get(hObject, 'String'));
if isnan(flip_thresh)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
% Hide (toggle callback) table if check button was checked or unchecked
if  get(handles.togglebutton1,'Value')== 1
    set(handles.togglebutton1,'Value',0); 
    togglebutton1_Callback(handles.togglebutton1, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function flip_thresh_CreateFcn(hObject, eventdata, handles)


% --- Export selected roi in workspace
function pushbutton_export_Callback(hObject, eventdata, handles)

if get(handles.togglebutton1,'Value')==0
    errordlg('Use show table button to select a region to explore','Error');
else
    % Gets selected items in java table
    all_idx = handles.table.jtable.getSelectedRows ; 
    idx = [] ;
    % Get proper indices in case table is sorted
    for kk=1:length(all_idx); idx(kk) = str2num(handles.table.jtable.getValueAt(all_idx(kk),handles.table.jtable.getColumnCount - 1)) ; end
    if isempty(idx)
        errordlg('You must select a region to export','Error');
    % Several items selected in the tablex
    elseif length(idx)>0
        fprintf(sprintf('Export %d regions(s) as variable "roi"\n',length(idx)));
        assignin('base','roi',handles.rois(idx))
   end
end


% --- Executes on button press in checkbox_significant_only.
function checkbox_significant_only_Callback(hObject, eventdata, handles)
% Hide (toggle callback) table if check button was checked or unchecked
if  get(handles.togglebutton1,'Value')== 1
    set(handles.togglebutton1,'Value',0); 
    togglebutton1_Callback(handles.togglebutton1, eventdata, handles);
end
