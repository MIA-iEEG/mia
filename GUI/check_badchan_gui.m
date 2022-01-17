function varargout = check_badchan_gui(varargin)
% CHECK_BADCHAN_GUI MATLAB code for check_badchan_gui.fig
%      CHECK_BADCHAN_GUI, by itself, creates a new CHECK_BADCHAN_GUI or raises the existing
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
    'gui_OpeningFcn', @check_badchan_gui_OpeningFcn, ...
    'gui_OutputFcn',  @check_badchan_gui_OutputFcn, ...
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


% --- Executes just before check_badchan_gui is made visible.
function check_badchan_gui_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for check_badchan_gui
handles.output = hObject;

% Loads good channels list from previous execution
goodlabels = varargin{1};
handles.init_labels=goodlabels;
set(handles.good_chan,'String',goodlabels);

% Loads bad channels list from previous execution
badlabels = varargin{2};
set(handles.bad_chan,'string',badlabels);

% Allows multiple selection in the list
set(handles.bad_chan,'Max',length(goodlabels));
set(handles.good_chan,'Max',length(goodlabels));

handles.Ogoodlabels=goodlabels ;
handles.Obadlabels=badlabels;

set(handles.figure1,'CloseRequestFcn',@my_closereq);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sanity_check_gui wait for user response (see UIRESUME)
%permet de faire attendre la sortie de ce programme ?? sanity_check
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = check_badchan_gui_OutputFcn(hObject, eventdata, handles)

% Cancel if user close window 
if isempty(handles)
    
    varargout{1} = {};
    varargout{2} = {};
    
else
    varargout{1} =cellstr(get(handles.good_chan,'string'));
    varargout{2} =cellstr(get(handles.bad_chan,'string'));
    
    % Close window 
    delete(handles.figure1);
end

% --- Executes during object creation, after setting all properties.
function good_chan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to good_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function bad_chan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bad_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in good_chan.
function good_chan_Callback(hObject, eventdata, handles)

% --- Executes on selection change in bad_chan.
function bad_chan_Callback(hObject, eventdata, handles)

% --- Executes on button press in BAD.
function BAD_Callback(hObject, eventdata, handles)
% hObject    handle to BAD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%on charge le nom des labels + on stock les listes initiales pour pouvoir
%les utiliser ?? l'appel du bouton cancel

Opening_goodlabels=get(handles.good_chan,'string');

Vec=zeros(1,length(Opening_goodlabels));
good2bad=get(handles.good_chan,'value');

%Echange des noms entre les deux listes
for ii=1:length(good2bad)
    Vec=Vec+strcmp(Opening_goodlabels,Opening_goodlabels(good2bad(ii)))';
end

set(handles.good_chan,'Value',1);

%actualisation de la liste des good et bad labels
goodlabels=Opening_goodlabels(find(~Vec));

%badlabels=sort(cat(1,get(handles.bad_chan,'string'),labels(find(Vec))));
badlabels=cat(1,get(handles.bad_chan,'string'),Opening_goodlabels(find(Vec)));

%stockage des donn??es
set(handles.good_chan,'string',goodlabels);
set(handles.bad_chan,'string',badlabels);


guidata(hObject,handles)



% --- Executes on button press in GOOD.
function GOOD_Callback(hObject, eventdata, handles)
% hObject    handle to GOOD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%on charge le nom des labels + on stock les listes initiales pour pouvoir
%les utiliser ?? l'appel du bouton cancel

Opening_badlabels=get(handles.bad_chan,'string');

Vec=zeros(1,length(Opening_badlabels));
bad2good=get(handles.bad_chan,'value');

for ii=1:length(bad2good)
    Vec=Vec+strcmp(Opening_badlabels,Opening_badlabels(bad2good(ii)))';
end

set(handles.bad_chan,'Value',1);

goodlabels=cat(1,get(handles.good_chan,'string'),Opening_badlabels(find(Vec)));
 % Reorder labels
        [~,res_index] = mia_sort_nat(goodlabels, 'ascend');
        goodlabels = goodlabels(res_index) ;
        badlabels=Opening_badlabels(find(~Vec));

set(handles.good_chan,'string',goodlabels);
set(handles.bad_chan,'string',badlabels);

guidata(hObject,handles)

% --- Executes on button press in OK.
function pushbutton_ok_Callback(hObject, eventdata, handles)
% hObject    handle to OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%fin de l'attente de sanity_check
uiresume(handles.figure1);

% --- Executes on button press in cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%on re-initialise les listes en cas d'annulation par l'utilisateur

set(handles.good_chan,'String',handles.Ogoodlabels);
set(handles.bad_chan,'string',handles.Obadlabels);

%fin de l'attente de sanity_check
uiresume(handles.figure1);

function my_closereq(src,callbackdata)
    errordlg('Use "cancel" button to close this figure');
