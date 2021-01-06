function varargout = sanity_check_gui(varargin)
% sanity_check_gui MATLAB code for sanity_check_gui.fig
%      SANITY_CHECK_GUI, by itself, creates a new SANITY_CHECK_GUI or raises the existing
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

gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sanity_check_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @sanity_check_gui_OutputFcn, ...
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

% --- Executes just before sanity_check_gui is made visible.
function sanity_check_gui_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for sanity_check_gui
handles.output = hObject;

handles.fname = cell2mat(varargin{1}) ; 

handles = initialize_gui(hObject, handles, false);

% Update handles structure
guidata(hObject, handles);


% --- Initialize the GUI
function handles = initialize_gui(fig_handle, handles, isreset)

% New File structure : optimization to load only average 
if ismember('Favg',who('-file', handles.fname))&&ismember('df',who('-file', handles.fname))
    % Loads data to display
    d = load(handles.fname,'Time','labels','Favg','df') ;
  
else
    % FIX OLD FILES
    d = load(handles.fname,'Time','labels','F') ;
    d.Favg = mean(d.F,3); Favg = d.Favg ;
    d.df = size(d.F,3) ; df = d.df ;
    save(handles.fname,'-append','Favg','df');
end

% Get patient name
[p,~,~] = fileparts(handles.fname)   ;
[outdir,pt,~] = fileparts(p)   ;

handles.Pt_name = pt;
handles.Fmono = d.Favg ;
handles.Labelsmono= d.labels ;
handles.df =d.df ; 

% Compute once for all bipolar montage
[handles.Fbi, handles.Labelsbi, handles.idx1, handles.idx2] = mia_make_bipolarmtg(d.Favg,d.labels)   ;

handles.Time = d.Time;
handles.outdir = outdir ;

handles.FONTSZ = 8 ; 

% Loading old file if they exist 
if ismember('isGood',who('-file', handles.fname))
    load(handles.fname,'isGood') ;
    handles.isNew = 0 ;
    handles.iSelMono =logical(isGood) ; 
    bad_chan = find(isGood==0) ; 
    iSel1 = ~sum(repmat(handles.idx1,length(bad_chan),1)==repmat(bad_chan,length(handles.idx1),1)',1) ; 
    iSel2 = ~sum(repmat(handles.idx2,length(bad_chan),1)==repmat(bad_chan,length(handles.idx2),1)',1) ; 
    handles.iSelBi = iSel1&iSel2;
           
else 
    handles.isNew = 1 ;
    handles.iSelMono = true(1,size(handles.Fmono,1)) ; 
    handles.iSelBi = true(1,size(handles.Fbi,1)) ; 

end 

handles.monobi=1;
handles.monobi_max=1;

update_butterfly(handles) ; 
update_imagesc(handles);
update_infos(handles) ;

% --- Outputs from this function are returned to the command line.
function varargout = sanity_check_gui_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
  

% --- Executes during object creation, after setting all properties.
function axes_imagesc_CreateFcn(hObject, eventdata, handles)

% --- Executes on selection change in popupmenu_montage.
function popupmenu_montage_Callback(hObject, eventdata, handles)

monobi=get(hObject,'value');
monobi_max=get(hObject,'Max');

if (monobi ~= monobi_max); 
    set(handles.SaveChannelsSelection,'enable','off'); 
%     handles.iSel = logical(ones(1,size(handles.Fbi,1))); 
else
    set(handles.SaveChannelsSelection,'enable','on'); 
%     handles.iSel = logical(ones(1,size(handles.Fmono,1))); 
end
     
update_butterfly(handles); 
update_imagesc(handles); 
update_infos(handles) ;

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_montage_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'Monopolar','Bipolar'});


% --- Executes during object creation, after setting all properties.
function uipanel_infos_CreateFcn(hObject, eventdata, handles)

% --- Get signal in proper montage
function [F, labels, iSel] = get_signals(handles) 

% Get Montage
if (get(handles.popupmenu_montage,'Value') == get(handles.popupmenu_montage,'Max'))
    F = handles.Fmono; 
    labels = handles.Labelsmono; 
    iSel = handles.iSelMono; 
    
else 
    F = handles.Fbi; 
    labels = handles.Labelsbi; 
    iSel = handles.iSelBi; 
end

% --- Update display of the butterfly signal view
function update_butterfly(handles) 

cla(handles.axes_average_butterfly);

% Get signal at the appropriate montage
[F,labels,iSel] = get_signals(handles) ;

% Get the good channels (as selected by user)
F = F(iSel,:,:);
labels = labels(iSel);

% Display TIMESERIES avergae : use this trick for older Matlab version (wont
% handle heterogeneous colormap on a same figure
clr=hsv(numel(unique(labels))); %get different color for each electrodes
Fmean=mean(F,3)';

for iVec=1:length(labels)
    
    plot(handles.axes_average_butterfly,handles.Time,Fmean(:,iVec),'color',clr(iVec,:)); %plot the signal of electrods
    hold(handles.axes_average_butterfly,'on');

end

grid(handles.axes_average_butterfly,'on') ; 
xlim(handles.axes_average_butterfly,[handles.Time(1),handles.Time(end)]);

% Set amplitude scale
m = max(max(abs(mean(F,3)))) ; 

ylim(handles.axes_average_butterfly,[-m,m]);

% Set Fontsize
set(handles.axes_average_butterfly, 'Fontsize',handles.FONTSZ);

xlabel(handles.axes_average_butterfly,'Time (ms)');
ylabel(handles.axes_average_butterfly,'Amplitude (uV)');

% Update signals on image (rigth axes)
function update_imagesc(handles) 

cla(handles.axes_imagesc);

% Get signal at the appropriate montage
[F,labels,iSel] = get_signals(handles) ;

% Get the good channels (as selected by user)
F = F(iSel,:,:);
labels = labels(iSel);

% Display IMAGE avergae
imagesc(handles.Time,1:size(F,1),mean(F,3), 'Parent',handles.axes_imagesc);

% Set amplitude scale
m = max(max(abs(mean(F,3)))) ; 
caxis(handles.axes_imagesc,[-m m]) ; 

colorbar('peer',handles.axes_imagesc,'location', 'NorthOutside');
colormap(handles.axes_imagesc,jet);

grid(handles.axes_imagesc);
set(handles.axes_imagesc,...
    'YTick',1:size(F,1),...
    'YTickLabel', strrep(labels,'_','\_'),...
    'Fontsize',handles.FONTSZ);
xlabel(handles.axes_imagesc,'Time (ms)');

% Update text infos
function update_infos(handles) 

% Get signal at the appropriate montage
[F,labels] = get_signals(handles) ;

% Fill out text fields 
set(handles.text_pt_name,'String',handles.Pt_name); 

% Display the number of trials ( = degree of freedom) 
set(handles.text_nb_trials,'String',num2str(handles.df)); 
pop_montage = get(handles.popupmenu_montage,'String') ; 

% Display the number of contacts 
set(handles.text_nb_contacts,'String', sprintf('%s (%s)',num2str(size(F,1)),cell2mat(pop_montage(get(handles.popupmenu_montage,'Value')))));

% Get labels of electrodes (NOT contacts)
elec = regexprep(labels,'\d+','') ; 
uelec = unique(elec);

set(handles.text_nb_electrodes,'String',num2str(length(uelec))); 

text_electrodes = '' ;

% Get the number of contact per electrodes and create list for diaply
for ii=1:length(uelec)
    nb_cont(ii) = sum(strcmp(uelec(ii),regexprep(labels,'\d+','')));
    text_electrodes = sprintf('%s\n%s \t:\t %d contacts',text_electrodes,uelec{ii},nb_cont(ii));
end

% Display list of electrodes with number of contacts 
set(handles.list_electrodes,'String', text_electrodes);


% --- Executes on button press in SelectChannels.
function SelectChannels_Callback(hObject, eventdata, handles)

% Get signal at the appropriate montage
[F,labels,iSel] = get_signals(handles) ;

[goodlabels,badlabels]=check_badchan_gui(labels(iSel),labels(~iSel));

% if monopolar montage 
if (get(handles.popupmenu_montage,'Value') == get(handles.popupmenu_montage,'Max'))

    % get 0 1 vector for labels that were marked as good
    handles.iSelMono=ismember(labels,goodlabels) ; 

    bad_chan = find(handles.iSelMono==0) ; 
    iSel1 = ~sum(repmat(handles.idx1,length(bad_chan),1)==repmat(bad_chan,length(handles.idx1),1)',1) ; 
    iSel2 = ~sum(repmat(handles.idx2,length(bad_chan),1)==repmat(bad_chan,length(handles.idx2),1)',1) ; 
    handles.iSelBi = iSel1&iSel2;
else
    
    handles.iSelBi=ismember(labels,goodlabels) ; 
end
update_butterfly(handles);
update_imagesc(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press SaveChannelsSelection Channel selection 
function SaveChannelsSelection_Callback(hObject, eventdata, handles)

% Construct a questdlg with three options
savedlg = questdlg('Are you sure you want to remove bad channels ?', ...
    'Save','Yes','No','Yes');

% Handle response
switch savedlg
    case 'Yes'
        if handles.isNew 
            % Appends isGood to the original file 
            isGood = handles.iSelMono ;
            save(handles.fname,'isGood','-append');
            handles.isNew =0;
        else
            choice = questdlg('This action will erase already existing file. Would you like to continue ? ', ...
                'CAUTION','Continue','Cancel','Continue');
            switch choice
                case 'Continue' 
                    isGood = handles.iSelMono ;
                    save(handles.fname,'isGood','-append');
                case 'Cancel'
                    
            end
        end
    case 'No'        
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press SaveChannelsSelection JPEG
function jpeg_export_ClickedCallback(hObject, eventdata, handles)
set(handles.uitoggletool1,'State','off');
set(handles.uitoggletool2,'State','off');

% Get default output directory and filename
[PATH,~,~]=fileparts(handles.outdir);

snap_filename=char(fullfile(PATH,'JPEGs'));

if ~exist(snap_filename,'dir')
    mkdir(snap_filename);
end

snap_filename=char(fullfile(snap_filename,strcat(handles.Pt_name,'_Sanity_Check.jpg')));

% Browse file to savechannelsselection
[filename,snap_filename]=uiputfile({'*.jpg;','Image Files';...
          '*.*','All Files' },'Save Image',...
          snap_filename);

% If savechannelsselection operation is cancelled
if snap_filename(1) == 0 ;  return; end 
    
snap_filename=char(fullfile(snap_filename,filename));     
       
% SaveChannelsSelection jpeg
export_fig(snap_filename,'-jpeg',handles.figure1)

guidata(hObject, handles);

% Zoom in Callback (do nothing, guide handles the action)
function uitoggletool1_OnCallback(hObject, eventdata, handles)

% Zoom out call back (do nothing, guide handles the action)
function uitoggletool2_OnCallback(hObject, eventdata, handles)
