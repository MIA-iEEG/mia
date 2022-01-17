function varargout = create_a_study_gui(varargin)
% GROUP_ANALYSIS_GUI MATLAB code for group_analysis_gui.fig
%      GROUP_ANALYSIS_GUI, by itself, creates a new GROUP_ANALYSIS_GUI or raises the existing
%      singleton*.
%
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
% 
% 2016/10/14 : ASD fixing the indexing issue in the java table (when
% filters are used)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @create_a_study_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @create_a_study_gui_OutputFcn, ...
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

% --- Executes just before group_analysis_gui is made visible.
function create_a_study_gui_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for statistics_gui
handles.output = hObject;
handles.idx_selected = varargin{1};
handles.outdir = varargin{2};
handles.current_loctable = varargin{3};
handles.INDEX = 8 ;

% Move window to the center of the screen 
movegui(gcf,'center');

%  ASD : create a copy of the table  (could be optimized??)
% Read the working directory in order to build the table
[handles.table.mia_table,handles.table.sFiles] = create_table_workdir(handles.outdir, handles.current_loctable) ;

jtable = com.jidesoft.grid.SortableTable(handles.table.mia_table,{'Patient','Method','Montage','Freq. band','Fs','Remove Avg','Nb stats','Localized Contacts','ID'});

% Trick to hide to indexing column
jtable.getColumnModel().getColumn(handles.INDEX).setMinWidth(0);
jtable.getColumnModel().getColumn(handles.INDEX).setMaxWidth(0);
jtable.getColumnModel().getColumn(handles.INDEX).setWidth(0);

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
jScrollPane = javax.swing.JScrollPane(jtable);

handles.table.jtable = jtable ;
handles.table.jScrollPane = jScrollPane ;
[handles.table.hjtable,handles.table.hjcontainer]=javacomponent(jScrollPane,[15,65,handles.figure1.Position(3)-20,handles.figure1.Position(4)-80],handles.figure1);

% Get intervals of consectuive rows selected in the table
a = handles.idx_selected;
deb= [a(1);a(find(diff(a)~=1)+1)] ;
fin = [a(find(diff(a)~=1));a(end)] ;

% It only works by giving intervals of selection..
for ss=1:length(deb)
    jtable.addRowSelectionInterval(deb(ss),fin(ss));
end

% data initialization
set(handles.nboot,'string','500');
set(handles.alpha,'string','0.001');
set(handles.low_twin,'string','-0.500');
set(handles.up_twin,'string','1.5');
set(handles.baselineb,'string','-0.8');
set(handles.baselinee,'string','-0.05' );

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes group_analysis_gui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = create_a_study_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1} = [] ;

% % The GUI is no longer waiting, just close it
% delete(handles.figure1);

% % --- Executes when user attempts to close figure1.
% function figure1_CloseRequestFcn(hObject, eventdata, handles)
% 
% hFig = ancestor(hObject,'Figure');
% if isequal(get(hFig,'waitstatus'),'waiting')
%     uiresume(hFig);
% else
%     delete(hFig);
% end


function nboot_Callback(hObject, eventdata, handles)
% hObject    handle to nboot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nboot as text
% str2double(get(hObject,'String')) returns contents of nboot as a double
nboot = str2double(get(hObject, 'String'));
if isnan(nboot)
    set(hObject, 'String', 500);
    errordlg('Input must be a number','Error');
elseif nboot<=0
    set(hObject, 'String', 500);
    errordlg('Input must be a positive value ','Error');
    
end

% Save the new nboot value
% handles.nboot = str2num(get(hObject,'string'));
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function nboot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nboot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha_Callback(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha as text
%str2double(get(hObject,'String')) returns contents of alpha as a double

alpha = str2double(get(hObject, 'String'));
if isnan(alpha)
    set(hObject, 'String', 0.001);
    errordlg('Input must be a number','Error');
elseif alpha<0 || alpha>1
    set(hObject, 'String', 0.001);
    errordlg('Input must be a value between 0 and 1','Error');
end

% Save the new alpha value
% handles.alpha = str2num(get(hObject,'string'));
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function low_twin_Callback(hObject, eventdata, handles)
% hObject    handle to low_twin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of low_twin as text
%str2double(get(hObject,'String')) returns contents of low_twin as a double
low_twin = str2double(get(hObject, 'String'));
if isnan(low_twin)
    set(hObject, 'String', -0.200);
    errordlg('Input must be a number','Error');
    
end

% Save the new alpha value
% handles.low_twin = str2num(get(hObject,'string'));
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function low_twin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to low_twin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function up_twin_Callback(hObject, eventdata, handles)
% hObject    handle to up_twin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of up_twin as text
%        str2double(get(hObject,'String')) returns contents of up_twin as a double
up_twin = str2double(get(hObject, 'String'));
if isnan(up_twin)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new alpha value
% handles.up_twin = str2num(get(hObject,'string'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function up_twin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to up_twin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calculate.
function calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Gets selected items in java table
all_idx = handles.table.jtable.getSelectedRows ; 

% Get proper indices in case table is sorted
for kk=1:length(all_idx); idx(kk) = str2num(handles.table.jtable.getValueAt(all_idx(kk),handles.INDEX)) ; end

mia_table=handles.table.mia_table(idx,:,:,:,:);

%check if selected lines have the same parameters method, freqid, montage
[method,~,~]=unique(mia_table(:,2));%list of all methods of selected files
[montage,~,~]=unique(mia_table(:,3));%list of all montage of selected fiels
[freq,~,~]=unique(mia_table(:,4));%list of all freq band of selected files
[srate,~,~]=unique(mia_table(:,5));%list of all freq band of selected files
[rmEvk,~,~]=unique(mia_table(:,6));%list of all freq band of selected files

%check that there is only 1 method, Motage and Freq.band through all
%selected patients
if ~(length(method)==1 && length(montage)==1 && length(freq)==1 && length(srate)==1 && length(rmEvk)==1)
    errordlg('You must select files with same parameters method, montage and freq. band','Files Selection');
    return;
end

% Folder GroupAnalysis creation if it doesn't exist yet
[PATH,NAME,EXT]=fileparts(handles.outdir);
outdir=fullfile(PATH,'GA_Results');

if ~exist(outdir,'dir')
    mkdir(outdir);
end

%%Ask for a file name to save data and verification if the given name alerady
%exist or not

%get files name already existing in outdir
d = dir(outdir);
isub = [d(:).isdir]; 
NAMES={d(isub).name}';
NAMES(ismember(NAMES,{'.','..'})) = [];

if isempty(NAMES)
    files=[];
else
    %get list of ganalysis files already existing 
    files={};
   for ii=1:length(NAMES)
       d=dir(fullfile(outdir,NAMES{ii}));      
       files=cat(1,files,{d.name}');
   end
  files(ismember(files,{'.','..'})) = []; 
  files(~logical(cellfun(@isempty,strfind(files,'cfroi'))))=[];

end


%function to check if the given name already exist
[NAME]=check_name(NAMES,outdir);
if isempty(NAME)
    return;
end

ganalysisOPTIONS.nboot=str2num(get(handles.nboot,'String'));
ganalysisOPTIONS.alpha=str2num(get(handles.alpha,'String'));
ganalysisOPTIONS.twin=[str2num(get(handles.low_twin,'String')), str2num(get(handles.up_twin,'String'))];
ganalysisOPTIONS.outdir=handles.outdir;
ganalysisOPTIONS.baseline=[str2num(get(handles.baselineb,'String')), str2num(get(handles.baselinee,'String'))] ;
ganalysisOPTIONS.subjects = handles.table.mia_table(idx,1) ; 

% Create study dir
mkdir(char(fullfile(outdir,NAME)));

%generate and save ganalysis = table of structure containing time and
%statistical data about selected subjects for the current study
ganalysis = s5_group_data_v1(handles.table.sFiles(idx),ganalysisOPTIONS);

% ASD 2018/1/26 : comment this line (at that point we don<t need the
% localisation table)
% [handles.rois,handles.table.mia_table,handles.edges,handles.datafiles] = create_struct_of_rois(handles.m_table_as, ganalysis, OPTIONS) ;

%save ganalysis
[subj,~,~]=unique(mia_table(:,1));

fname =char(fullfile(outdir,NAME,strcat(NAME,'_',num2str(length(subj)),'_',method,'_',montage,'_',freq)));
save(fname,'ganalysis');


% %save ganalysis
% [subj,~,~]=unique(mia_table(:,1));
% % Gets the frequency band for this processing 
% freqb = ganalysis{1}.freqb ;
% fname =char(fullfile(outdir,NAME,strcat(NAME,'_',num2str(length(subj)),'pts_',method,'_',montage,'_',strcat(num2str(freqb(1)),'_',num2str(freqb(end))))));
% save(fname,'ganalysis');



% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(handles.figure1);

function [studyn]=check_name(NAMES,outdir)

% Prompt user for a study name
prompt = {'Enter study name:'};
dlg_title = 'STUDY NAME';
num_lines = 1;
def = {'Study','hsv'};
studyn = mia_newid(prompt,dlg_title,num_lines,def);

%Use of cancel button
if isempty(studyn)
    return;
end

%files existing in outdir
if ~isempty(NAMES)
    %if the name already exist :
    if sum(ismember(NAMES,studyn))==1
        
        %creation of a list containing all the names that have been used in the
        %past
        list='';
        for ii=1:length(NAMES)
            list=[list, '  ',char(NAMES(ii))];
        end
        
        % ask for overwriting or choose an other one
        choice = questdlg(sprintf('Those names are alerady existing :%s',list) , ...
            'STUDY NAME','Overwrite','Choose an other Name','Overwrite');
        switch choice
            case 'Overwrite'
                %find back the full existing file name and delete it                
                 rmdir(char(fullfile(outdir,studyn)),'s');                   
                
            case 'Choose an other Name'
                %function to check if the given name already exist
                studyn=check_name(NAMES,outdir);
        end
    
    end    
end


function baselineb_Callback(hObject, eventdata, handles)
% hObject    handle to baselineb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of baselineb as text
%        str2double(get(hObject,'String')) returns contents of baselineb as a double
baselineb = str2double(get(hObject, 'String'));
if isnan(baselineb)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new alpha value
% handles.baselineb = str2num(get(hObject,'string'));
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function baselineb_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function baselinee_Callback(hObject, eventdata, handles)
baselinee = str2double(get(hObject, 'String'));
if isnan(baselinee)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

% Save the new alpha value
% handles.baselinee = str2num(get(hObject,'string'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function baselinee_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in allow_flipsign.
function allow_flipsign_Callback(hObject, eventdata, handles)
% hObject    handle to allow_flipsign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of allow_flipsign
