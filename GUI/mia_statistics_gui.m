function varargout = mia_statistics_gui(varargin)
% MIA_STATISTICS_GUI MATLAB code for mia_statistics_gui.fig
%      MIA_STATISTICS_GUI, by itself, creates a new MIA_STATISTICS_GUI or raises the existing
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
% 2016/10/14 : ASD fixing the indexing issue in the java table (when
% filters are used)

gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mia_statistics_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @mia_statistics_gui_OutputFcn, ...
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


% --- Executes just before mia_statistics_gui is made visible.
function mia_statistics_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mia_statistics_gui (see VARARGIN)

% Choose default command line output for mia_statistics_gui
handles.output = hObject;
idx_selected = varargin{1};
handles.outdir = varargin{2};
handles.current_loctable = varargin{3};
handles.INDEX = 8 ;

% Move window to the center of the screen 
movegui(gcf,'center');

%  ASD : create a copy of the table  (could be optimized??)
% Read the working directory in order to build the table
[handles.table.mia_table,handles.table.sFiles] = mia_create_table_workdir(handles.outdir, handles.current_loctable) ;

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
[handles.table.hjtable,handles.table.hjcontainer]=javacomponent(jScrollPane,[15,55,handles.figure1.Position(3)-20,handles.figure1.Position(4)-80],handles.figure1);

% % Get intervals of consectuive rows selected in the table
% a = idx_selected;
% deb= [a(1);a(find(diff(a)~=1)+1)] ;
% fin = [a(find(diff(a)~=1));a(end)] ;
% 
% % It only works by giving intervals of selection..
% for ss=1:length(deb) 
%     jtable.addRowSelectionInterval(deb(ss),fin(ss));
% end

% % Set edit BASELINE fields
% baseline_b = min(min(s.Time(s.Time<0))) ; 
% baseline_e = max(s.Time(s.Time<0)) ; 
% 
% set(handles.edit_baselineb,'String',num2str(baseline_b));
% set(handles.edit_baselinee,'String',num2str(baseline_e));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mia_statistics_gui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mia_statistics_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles) 
    varargout{1} = [] ; 
else
    
    % Get default command line output from handles structure
    varargout{1} = handles.output;
end


function edit_nboot_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nboot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nboot = str2double(get(hObject, 'String'));
if isnan(nboot)
    set(hObject, 'String', 500);
    errordlg('Input must be a number','Error');
elseif nboot<=0
    set(hObject, 'String', 500);
    errordlg('Input must be a positive value ','Error');

end

% Save the new up_freq value
handles.metricdata.nboot = nboot;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edit_nboot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nboot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to edit_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

alpha = str2double(get(hObject, 'String'));
if isnan(alpha)
    set(hObject, 'String', 0.001);
    errordlg('Input must be a number','Error');
elseif alpha<0 || alpha>1
    set(hObject, 'String', 0.001);
    errordlg('Input must be a value between 0 and 1','Error');
end

% Save the new up_freq value
handles.metricdata.alpha = alpha;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edit_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_calculate.
function pushbutton_calculate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idx = [] ; 


% Gets selected items in java table (+1 because java starts indexing at 0)
all_idx = handles.table.jtable.getSelectedRows ; 

% Get proper indices in case table is sorted
for kk=1:length(all_idx); idx(kk) = str2num(handles.table.jtable.getValueAt(all_idx(kk),handles.INDEX)) ; end

if isempty(idx)
    
    errordlg('You must select a stats file','Error');
   
else
    
    statOPTIONS.outdir= handles.outdir ;
    statOPTIONS.nboot = str2num(get(handles.edit_nboot,'String'));
    statOPTIONS.alpha =str2num(get(handles.edit_alpha,'String'));
    statOPTIONS.baseline =[str2num(get(handles.edit3,'String')),str2num(get(handles.edit4,'String'))] ;

    for ii=1:length(idx)
        mia_s5_compute_stats(handles.table.sFiles{idx(ii)},statOPTIONS);
    end
    
    files_list = sprintf('''%s'',',handles.table.sFiles{idx});
    
     % Stack in history
    mia_cmd_history(sprintf('MIA command : mia_s5_compute_stats(files,OPTIONS). \nfiles = \n%s\nOPTIONS = \n\toutdir = %s\n\tnboot= %s\n\talpha = %0.3f\n\tbaseline = [%0.3f,%0.3f]\n',...
                files_list,...
                statOPTIONS.outdir,...
                num2str(statOPTIONS.nboot),...
                statOPTIONS.alpha,...
                statOPTIONS.baseline(1),statOPTIONS.baseline(2)));

    close(handles.figure1);
end


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figure1);



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
baselineb = str2double(get(hObject, 'String'));
if isnan(baselineb)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
baselinee = str2double(get(hObject, 'String'));
if isnan(baselinee)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
