function varargout = extract_frequency(varargin)
% EXTRACT_FREQUENCY MATLAB code for extract_frequency.fig
%      EXTRACT_FREQUENCY, by itself, creates a new EXTRACT_FREQUENCY or raises the existing
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
                   'gui_OpeningFcn', @extract_frequency_OpeningFcn, ...
                   'gui_OutputFcn',  @extract_frequency_OutputFcn, ...
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

% --- Executes just before extract_frequency is made visible.
function extract_frequency_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for extract_frequency
handles.output = hObject;

%give a name to the figure
set(handles.figure1,'Name','Extract Frequency');

% Move window to the center of the screen 
movegui(gcf,'center');

handles.DEFAULTNCYCLES = 7;

list_patients = get(varargin{1},'String');
selected_patient = get(varargin{1},'Value'); 

if isempty(selected_patient)
    selected_patient  = 1 ; 
end
% Get list content and items selected
set(handles.listbox_patients,'Max',length(list_patients));
set(handles.listbox_patients,'Value',selected_patient);
set(handles.listbox_patients,'String',list_patients);

% Get working directory
handles.outdir = varargin{2} ;

% Loads first patient selected in list to set baseline default
tmp = dir(fullfile(handles.outdir,cell2mat(list_patients(selected_patient(1))),'*signal_LFP.mat')) ;      
s = load(cell2mat(fullfile(handles.outdir,list_patients(selected_patient(1)),tmp.name)),'Time');

% Set edit BASELINE fields
baseline_b = min(min(s.Time(s.Time<0))) ; 
baseline_e = max(s.Time(s.Time<0)) ; 

set(handles.edit_baselineb,'String',num2str(baseline_b));
set(handles.edit_baselinee,'String',num2str(baseline_e));

update_nPatient(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes extract_frequency wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = extract_frequency_OutputFcn(hObject, eventdata, handles)

% This line handles the case where user close the window 
if isfield(handles,'sFiles') ; varargout{1} = handles.sFiles ; else varargout{1} = [] ; end

% The GUI is no longer waiting, just close it
delete(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)

hFig = ancestor(hObject,'Figure');
if isequal(get(hFig,'waitstatus'),'waiting')
    uiresume(hFig);
else
    delete(hFig);
end

% --- Executes during object creation, after setting all properties.
function low_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to low_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function low_freq_Callback(hObject, eventdata, handles)
% hObject    handle to low_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of low_freq as text
%        str2double(get(hObject,'String')) returns contents of low_freq as a double
lfreq = str2double(get(hObject, 'String'));
if isnan(lfreq)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

update_edgeEffect(handles);

% Save the new low_freq value
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function up_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to up_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function up_freq_Callback(hObject, eventdata, handles)
% hObject    handle to up_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of up_freq as text
%        str2double(get(hObject,'String')) returns contents of up_freq as a double
ufreq = str2double(get(hObject, 'String'));
if isnan(ufreq)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new up_freq value
guidata(hObject,handles)

% --- Executes on button press in calculate.
function handles = calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Get montage
if (get(handles.radio_monopolar,'Value')) ; mtg = 'monopolar'; else ; mtg = 'bipolar'; end

% Get list of selected patients 
list_patients = get(handles.listbox_patients,'String');
selected_patient = get(handles.listbox_patients,'Value'); 

%% LFP extraction
contents = cellstr(get(handles.popupmenu_methods,'String')) ;

% Get montage, dir, subjects to process
OPTIONS.mtg = mtg;
OPTIONS.outdir=handles.outdir ; 
OPTIONS.subjects = list_patients(selected_patient) ;

% Get baseline bounds 
baseline_b = str2num(get(handles.edit_baselineb,'String')) ; 
baseline_e = str2num(get(handles.edit_baselinee,'String')) ;
OPTIONS.zbaseline = [baseline_b,baseline_e]; 
OPTIONS.overwrite = 'no'; % Overwrite any existing files

% if (get(handles.checkbox_overwrite,'Value')) ; OPTIONS.overwrite = 'yes'; else OPTIONS.overwrite = 'No'; end;
all_methods =get(handles.popupmenu_methods,'String')  ;
modetf = all_methods(get(handles.popupmenu_methods,'Value'));

if strcmpi(modetf,'Morse') %% COMPUTE LFP EXTRACTION
    
    %% Frequency extraction 
    % Test if frequency bounds are valid
    lfreq = str2num(get(handles.low_freq,'String')); 
    ufreq = str2num(get(handles.up_freq,'String')); 
    step =  str2num(get(handles.step,'String')); 
    
    if ufreq<=lfreq 
       warndlg('High Frequency bound must be > Low frequency bound') ;
       return ;
    end

    if isempty(step)||(step==0)
        warndlg('Missing value : Frequency step') ;
       return ;
    end

    %Options = Select descomposition mode, frequency band and baseline
    OPTIONS.modetf=modetf ; 
    OPTIONS.freqs= [lfreq:step:ufreq];
    OPTIONS.subjects = list_patients(selected_patient) ; 

    % Compute frequecy analysis
    handles.sFiles = mia_compute_morse(OPTIONS);
   
    
elseif strcmp(modetf,'LFP') %% COMPUTE LFP EXTRACTION
    
   handles.sFiles = mia_s3_zscore(OPTIONS) ; 
 
else %% COMPUTE TIME-FREQUENCY EXTRACTION
    
    %% Frequency extraction 
    % Test if frequency bounds are valid
    lfreq = str2num(get(handles.low_freq,'String')); 
    ufreq = str2num(get(handles.up_freq,'String')); 
    step =  str2num(get(handles.step,'String')); 
    removeEvoked = get(handles.checkbox_remove_evoked,'Value');

    if ufreq<=lfreq 
       warndlg('High Frequency bound must be > Low frequency bound') ;
       return ;
    end

    if isempty(step)||(step==0)
        warndlg('Missing value : Frequency step') ;
       return ;
    end

    %Options = Select descomposition mode, frequency band and baseline
    OPTIONS.modetf=get(handles.popupmenu_methods,'String');

    OPTIONS.freqs= [lfreq:step:ufreq];
    OPTIONS.subjects = list_patients(selected_patient) ; 
    OPTIONS.removeEvoked = removeEvoked;
    OPTIONS.ncycles=str2num(get(handles.edit_ncycles,'String'));

    % Compute frequecy analysis
    handles.sFiles = mia_s4_compute_tf_bandwise(OPTIONS);
   
end

% Close GUI
guidata(hObject, handles);
figure1_CloseRequestFcn(hObject,handles);

% --- Executes on button press in button_cancel.
function button_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to button_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure1_CloseRequestFcn(hObject,handles);


% --- Executes when selected object changed in unitgroup.
function unitgroup_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu_methods.
function popupmenu_methods_Callback(hObject, eventdata, handles)

% Get current value in the popupmenu : Morlet? LFP?
contents = cellstr(get(hObject,'String')) ;
method =contents{get(hObject,'Value')}  ;

% If LFP : disable all frequency related edit fields
if strcmp(method,'LFP')
    set(handles.low_freq, 'String', '0');
    set(handles.up_freq,  'String', '0');
    set(handles.step,     'String', '0');     
    set(handles.edit_ncycles,     'String', '0');     
    set(handles.text_edges,'String','');

    set(handles.low_freq, 'enable', 'off');
    set(handles.up_freq,  'enable', 'off');
    set(handles.step,     'enable', 'off');     
    set(handles.edit_ncycles,     'enable', 'off');     
    set(handles.checkbox_remove_evoked, 'enable','off');
else
    set(handles.low_freq, 'enable', 'on');
    set(handles.up_freq,  'enable', 'on');
    set(handles.step,     'enable', 'on');     
    set(handles.edit_ncycles,     'enable', 'on');
    set(handles.checkbox_remove_evoked, 'enable','on');
end
%  contents{get(hObject,'Value')} returns selected item from popupmenu_methods


% --- Executes during object creation, after setting all properties.
function popupmenu_methods_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_methods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% set(hObject,'String',{'Morlet wavelet';'LFP';'Morse'});
set(hObject,'String',{'Morlet wavelet';'LFP'});

% --- Executes on selection change in listbox_patients.
function listbox_patients_Callback(hObject, eventdata, handles)

update_edgeEffect(handles) ;     
update_nPatient(handles) ; 
guidata(hObject,handles) ;


% --- Update NPatient text filed as function of selected patietns in list
function update_nPatient(handles)
        
selected_patient = get(handles.listbox_patients,'Value'); 

set(handles.nPatients,'String',num2str(length(selected_patient)));


% --- Update NPatient text filed as function of selected patietns in list
function update_edgeEffect(handles)
   
%Get time information of selected patient
list_patients = get(handles.listbox_patients,'String');
selected_patient = get(handles.listbox_patients,'Value');

if ~isempty(selected_patient) & str2num(get(handles.low_freq,'String'))~= 0;
        
    tmp = dir(cell2mat(fullfile(handles.outdir,list_patients(selected_patient(1)),'*signal_LFP.mat'))) ;
    s = load(cell2mat(fullfile(handles.outdir,list_patients(selected_patient(1)) ,tmp.name)),'Time');
    
    lfreq=str2num(get(handles.low_freq,'String')); 
        
    ncycles = str2num(get(handles.edit_ncycles,'String'));
    
    % Compute baseline (calculate edges effects)
    baseline_b=min(min(s.Time(s.Time<0)))+ncycles./(pi*lfreq);
    baseline_e = max(s.Time(s.Time<0)) - ncycles./(pi*lfreq) ; 
   
    set(handles.edit_baselineb,'String',num2str(baseline_b));
    set(handles.edit_baselinee,'String',num2str(baseline_e));

    %if wrong baseline (baseline_b>baseline_e) => display warning message
    if baseline_b>baseline_e
        warndlg('Baseline Contain Edges Effect. Play at your own risk');
        time=s.Time(s.Time<0);
        baseline_b=time(1);
        baseline_e = time(end); 
    end
    
    textebaseline=sprintf('(Edges effects : %0.3f sec)',ncycles./(pi*lfreq) );
    set(handles.text_edges,'String',textebaseline);

    
end


% --- Executes during object creation, after setting all properties.
function listbox_patients_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_patients (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_baselineb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_baselineb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_baselineb_Callback(hObject, eventdata, handles)
% hObject    handle to edit_baselineb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_baselineb as text
%        str2double(get(hObject,'String')) returns contents of edit_baselineb as a double
% --- Executes during object creation, after setting all properties.


% --- Executes during object creation, after setting all properties.
function edit_baselinee_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_baselinee (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_baselinee_Callback(hObject, eventdata, handles)
% hObject    handle to edit_baselineb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_baselineb as text
%        str2double(get(hObject,'String')) returns contents of edit_baselineb as a double
% --- Executes during object creation, after setting all properties.



function step_Callback(hObject, eventdata, handles)
% hObject    handle to step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of step as text
%        str2double(get(hObject,'String')) returns contents of step as a double
step = str2double(get(hObject, 'String'));
if isnan(step)
    set(hObject, 'String', 1);
    errordlg('Input must be a number','Error');
end

if step==0
     set(hObject, 'String', 1);
    errordlg('Invalid value for step','Error');
end

if step>(str2num(get(handles.up_freq,  'String'))-str2num(get(handles.low_freq,  'String')))
      set(hObject, 'String', 1);
    errordlg('Invalid value for step','Error');
end

% Save the new up_freq value
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to up_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_methods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_ncycles_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ncycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ncycles as text
%        str2double(get(hObject,'String')) returns contents of edit_ncycles as a double
ncylces = str2double(get(hObject, 'String'));
if isnan(ncylces)
    set(hObject, 'String', handles.DEFAULTNCYCLES);
    errordlg('Input must be a number','Error');
end

update_edgeEffect(handles) ;

% Save the new low_freq value
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function edit_ncycles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ncycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text_edges_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_edges (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function radio_bipolar_CreateFcn(hObject, eventdata, handles)

function radio_bipolar_Callback(hObject, eventdata, handles)

function radio_monopolar_Callback(hObject, eventdata, handles)


% --- Executes on button press in checkbox_remove_evoked.
function checkbox_remove_evoked_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_remove_evoked (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_remove_evoked
