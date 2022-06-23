function varargout = mia_display_images_stats_gui(varargin)
% MIA_DISPLAY_IMAGES_STATS_GUI MATLAB code for mia_display_images_stats_gui.fig
%      MIA_DISPLAY_IMAGES_STATS_GUI, by itself, creates a new MIA_DISPLAY_IMAGES_STATS_GUI or raises the existing
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
% Copyright (C) 2016-2022 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mia_display_images_stats_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @mia_display_images_stats_gui_OutputFcn, ...
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


% --- Executes just before mia_display_images_stats_gui is made visible.
function mia_display_images_stats_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mia_display_images_stats_gui (see VARARGIN)

% Choose default command line output for mia_display_images_stats_gui
handles.output = hObject;
 
% This sets up the initial plot - only do when we are invisible
handles.fname = cell2mat(varargin{1}) ;
handles.infos = varargin{2} ;

%initilization datas for select channels button 
handles.badlabels={};

% Move window to the center of the screen 
movegui(hObject,'center');
set(hObject,'color','white');

% Initialize GUI 
handles = initialize_gui(handles);

% Track mouse cursor in order to display latencies in real time
set (hObject, 'WindowButtonMotionFcn', {@mouseMove, handles});

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function handles = initialize_gui(handles)
% initialize the whole GUI 

handles.FONTSZ = 8 ; 
handles.FONTSZ_large = 14 ; 

handles.BACKGROUNDCOL = [ 0.8275, 0.8275, 0.8275] ; 

% New File structure : optimization to load only average 
if ismember('meanzs',who('-file', handles.fname))&&ismember('tvals',who('-file', handles.fname))&&ismember('pvals',who('-file', handles.fname))&&ismember('df',who('-file', handles.fname))
    % Loads data to display
    d = load(handles.fname,'Time','labels','meanzs','tvals','pvals','df') ;
  
else
    % FIX OLD FILES
    d = load(handles.fname,'Time','labels','zs') ;
    d.meanzs = mean(d.zs,3); meanzs = d.meanzs ;
    d.df = size(d.zs,3) ; df = d.df ;
    % Compute the t-test
    [d.tvals,d.pvals] = mia_compute_ttest(d.zs);
    tvals = d.tvals ; pvals = d.pvals;
    % Fix old files
    save(handles.fname,'-append','meanzs','tvals','pvals','df');
end

% Check if stats file exists 
[PATHSTR,NAME,EXT] = fileparts(handles.fname) ; 
fname_stat = fullfile(PATHSTR,strcat(strrep(NAME,'_data','_stats'),EXT)) ;

% handles.zs = zs ;
handles.Labels = d.labels;
handles.Time = d.Time ;
handles.tvals = d.tvals; 
handles.pvals = d.pvals; 
 
% ASD TODO : display name of regions if contacts were labelled (m_table
% exist)
if ~exist(fname_stat)
  hwarn=warndlg('No statsistics were computed for this data. Play at your own risk') ;
  waitfor(hwarn);
  
  handles.nboot = 0 ;
  handles.pthresh = 0.001 ;
  handles.threshdur = 0.02 ;
  
else
    stats = load(fname_stat) ;

    % Ensure compatibility with old files
    if isfield(stats,'stats')
          handles.nboot = [stats.stats(:).nboot] ;
          handles.pthresh = [stats.stats(:).pthresh];
          handles.threshdur = [stats.stats(:).threshdur] ;

   else
        % ASD TODO : test with old data
        Fs = 1/(Time(2)-Time(1));
        handles.nboot =-1 ;
        handles.pthresh = 0.001;
        handles.threshdur = stats.threshdur/Fs ;
    end

end

% zsnew used to handle channel selection
handles.meanzsnew=d.meanzs;

% Check if montage exist
fname_montage = fullfile(PATHSTR,strcat(NAME,'_montage.mat')) ; 
if exist(fname_montage,'file')
    load(fname_montage) ; 
    
    % Fix old files 
    if exist('isGood', 'var') 
        isDisplayed = isGood ; 
        save(fname_montage,'isDisplayed') ; 
    end
        
    handles.iSel = isDisplayed ;
else  
    handles.iSel = ones(1,size(handles.meanzsnew,1))==1 ; 
end

% Set default scale to maximum value
% Set the zscore editable text filed to [default max value]
max_zs = max(max(abs(handles.meanzsnew))) ;
set(handles.edit_scale,'String',sprintf('%0.2f',max_zs));

% Set the zscore static text to "Scale (default max value)"
set(handles.text_scalezs,'String',sprintf('Max z-score scale (%0.2f)',max_zs));

% Fill in the pop up menu with precompute stats
for ii=1:length(handles.pthresh)
   
    tmp{ii} = sprintf('Parameters : {nBoot=%d   p < %1.4f} Results : {duration = %1.4f}',handles.nboot(ii),handles.pthresh(ii),handles.threshdur(ii));
    
end

set(handles.popupmenu_preconfig,'String',tmp) ; 

% Init the editable fields for pval and duration 
set(handles.edit_pvalue,'String',num2str(handles.pthresh(1))) ; 
set(handles.edit_duration,'String',num2str(handles.threshdur(1))) ; 

% Set first line of text (top)
set(handles.text_file_infos,'String',sprintf('Patient : %s\t Method : %s\t Montage : %s\t  Freqs : %s\t',handles.infos{1},handles.infos{2},handles.infos{3},handles.infos{4}));

% Set the tvals editable text filed to [default max value]
max_tvals = max(max(abs(mean(d.tvals,3)))) ;
set(handles.edit_scale_tvals,'String',sprintf('%0.2f',max_tvals));

% Set the tvals static text to "Scale (default max value)"
set(handles.text_scaletv,'String',sprintf('Max t-value scale (%0.2f)',max_tvals));

% Display first column (left)
update_orig(handles);

% Display second column (middle)
update_stats(handles);

% Display third column (right)
update_stats_filtertime(handles);

% Link the three axes so that zoom is common to the three
ax(1)=handles.axes_duration;
ax(2)=handles.axes_stats;
ax(3)=handles.axes_orig;

linkaxes([ax(1) ax(2) ax(3)],'xy');

% Update handles structure
guidata(handles.figure1, handles);

% --- Display zscored values (power)
function varargout = update_orig(handles) 

% Get all parameters for visulalization (from edit fields)
[threshp, threshdur,max_tvals,max_zscore] = get_parameters(handles) ; 

time = handles.Time; 
meanzs = handles.meanzsnew(handles.iSel,:);

% Display IMAGE avergae

hImage = imagesc(time,1:size(meanzs,1),meanzs,'parent',handles.axes_orig);
caxis(handles.axes_orig,[-max_zscore max_zscore]) ; 
colorbar('peer',handles.axes_orig,'location', 'NorthOutside','Fontsize',handles.FONTSZ_large);
grid(handles.axes_orig);

% Display electrode labels if FONTSZ~=0 
if handles.FONTSZ~=0
    set(handles.axes_orig,...
    'YTick',1:size(meanzs,1),...
    'YTickLabel', strrep(handles.Labels(handles.iSel),'_','\_'),...
    'Fontsize',handles.FONTSZ);
else
    set(handles.axes_orig,'YTick',...
        1:size(meanzs,1),...
        'YTickLabel',''); 
end

xlabel(handles.axes_orig,'Time (sec)','Fontsize',handles.FONTSZ_large);

% Light grey color for middle of the scale 
cmap = colormap(jet) ; cmap(fix(length(cmap)/2)+1,:) = handles.BACKGROUNDCOL; colormap(handles.axes_orig,cmap) ;


% --- Display Statistic axes (tvals)
function varargout = update_stats(handles) 

% Get all parameters for visulalization (from edit fields)
[threshp, threshdur,max_tvals,max_zscore] = get_parameters(handles) ; 

tvals = handles.tvals(handles.iSel,:);
pvals = handles.pvals(handles.iSel,:);
time = handles.Time; 

% Compute mask 
h = pvals<threshp; 

% Display IMAGE avergae

hImage2 = imagesc(time,1:size(tvals,1),mean(tvals,3).*h, 'Parent',handles.axes_stats );
colorbar('peer',handles.axes_stats ,'location', 'NorthOutside','Fontsize',handles.FONTSZ_large);
 
% % Light grey color for middle of the scale 
% cmap = colormap ; cmap(fix(length(cmap)/2)+1,:) = handles.BACKGROUNDCOL; colormap(cmap) ;
cmap = colormap(jet) ; cmap(fix(length(cmap)/2)+1,:) = handles.BACKGROUNDCOL; colormap(handles.axes_stats,cmap) ;

caxis(handles.axes_stats,[-max_tvals   max_tvals]) ; 
grid(handles.axes_stats) ;

% Display electrode labels if FONTSZ~=0 
if handles.FONTSZ~=0
    set(handles.axes_stats,...
        'YTick',1:size(tvals,1),...
        'YTickLabel',strrep(handles.Labels(handles.iSel),'_','\_'),...
        'Fontsize',handles.FONTSZ);
else
    set(handles.axes_stats,'YTick',...
        1:size(tvals,1),...
        'YTickLabel',''); 
end

xlabel(handles.axes_stats,'Time (sec)','Fontsize',handles.FONTSZ_large);



% --- Display Statistic axes (tvals)
function varargout = update_stats_filtertime(handles) 

% Get all parameters for visulalization (from edit fields)
[threshp, threshdur,max_tvals,max_zscore] = get_parameters(handles) ; 

tvals = handles.tvals(handles.iSel,:);
pvals = handles.pvals(handles.iSel,:);
time = handles.Time; 

% Compute mask 
h = pvals<threshp; 

Fs = 1/(time(2)-time(1));
[hf] = mia_filter_timewin_signif(h,threshdur*Fs);
     
% Display IMAGE avergae
%modif by Jane
hImage2 = imagesc(time,1:size(tvals,1),mean(tvals,3).*hf,'Parent',handles.axes_duration);
colorbar('peer',handles.axes_duration ,'location', 'NorthOutside','Fontsize',handles.FONTSZ_large);
caxis(handles.axes_duration,[-max_tvals  max_tvals ]) ; 

% Light grey color for middle of the scale 
cmap = colormap(jet) ; cmap(fix(length(cmap)/2)+1,:) = handles.BACKGROUNDCOL; colormap(handles.axes_duration,cmap) ;

grid(handles.axes_duration) ;
 
% Display electrode labels if FONTSZ~=0 
if handles.FONTSZ~=0

    set(handles.axes_duration,...
        'YTick',1:size(tvals,1),...
        'YTickLabel',strrep(handles.Labels(handles.iSel),'_','\_'),...
        'Fontsize',handles.FONTSZ);
else
    set(handles.axes_duration,'YTick',...
        1:size(tvals,1),...
        'YTickLabel',''); 
end

xlabel(handles.axes_duration,'Time (sec)','Fontsize',handles.FONTSZ_large);



% --- Display Statistic axes (tvals)
function [threshp, threshdur,max_tvals,max_zscore] = get_parameters(handles) 

threshp = str2num(get(handles.edit_pvalue,'String')) ; 
threshdur = str2num(get(handles.edit_duration,'String')) ; 

max_tvals = str2num(get(handles.edit_scale_tvals,'String')) ; 
max_zscore = str2num(get(handles.edit_scale,'String')) ; 


% --- Outputs from this function are returned to the command line.
function varargout = mia_display_images_stats_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu_preconfig.
function popupmenu_preconfig_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_preconfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_preconfig contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_preconfig
idx =get(hObject,'Value') ; 

% Set the edit text box (p and duration)
set(handles.edit_pvalue,'String',num2str(handles.pthresh(idx)));
set(handles.edit_duration,'String',num2str(handles.threshdur(idx)));

% Update displays
update_stats_filtertime(handles);
update_stats(handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_preconfig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_preconfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pvalue_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_stats(handles);
update_stats_filtertime(handles);

% --- Executes during object creation, after setting all properties.
function edit_pvalue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_duration_Callback(hObject, eventdata, handles)
% hObject    handle to edit_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_duration as text
%        str2double(get(hObject,'String')) returns contents of edit_duration as a double

threshdur = str2double(get(hObject, 'String'));
if isnan(threshdur) | (length(threshdur)~=1)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
else

    % Update third column (right)
    update_stats_filtertime(handles);

end

% --- Executes during object creation, after setting all properties.
function edit_duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_scale_Callback(hObject, eventdata, handles)
% hObject    handle to edit_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_scale as text
%        str2double(get(hObject,'String')) returns contents of edit_scale as a double
val_max = str2double(get(hObject, 'String'));
if isnan(val_max)
    set(hObject, 'String', 0.001);
    errordlg('Input must be a number','Error');
else
    handles.val_max = val_max;
    caxis(handles.axes_orig,[-val_max val_max]) ; 
end

% --- Executes during object creation, after setting all properties.
function edit_scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axes_orig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_orig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




function edit_scale_tvals_Callback(hObject, eventdata, handles)
% hObject    handle to edit_scale_tvals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_scale_tvals as text
%        str2double(get(hObject,'String')) returns contents of edit_scale_tvals as a double
tval_max = str2double(get(hObject, 'String'));
if isnan(tval_max)
    set(hObject, 'String', 0 );
    errordlg('Input must be a number','Error');
else
    handles.tval_max = tval_max;
    caxis(handles.axes_stats,[-tval_max tval_max]) ; 
    caxis(handles.axes_duration,[-tval_max tval_max]) ; 
end

% --- Executes during object creation, after setting all properties.
function edit_scale_tvals_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_scale_tvals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in Select_chan.
function Select_chan_Callback(hObject, eventdata, handles)
% hObject    handle to Select_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[goodlabels,badlabels]=mia_check_badchan_gui(handles.Labels(handles.iSel),handles.Labels(~handles.iSel));

labels=handles.Labels;

% get 0 1 vector for labels that were marked as good
handles.isDisplayed=ismember(labels,goodlabels) ; 
handles.iSel =  handles.isDisplayed==1 ;

% Display first column (left)
update_orig(handles);

% Display second column (middle)
update_stats(handles);

% Display third column (right)
update_stats_filtertime(handles);

% Save Channel montage 
isDisplayed = handles.isDisplayed ;
[PATHSTR,NAME,EXT] = fileparts(handles.fname) ; 
fname = fullfile(PATHSTR,strcat(NAME,'_montage.mat')) ; 
save(fname,'isDisplayed'); 

% Save the new low_freq value
guidata(hObject,handles)


% --------------------------------------------------------------------
function jpeg_export_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to jpeg_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uitoggletool1,'State','off');
set(handles.uitoggletool2,'State','off');
set(handles.uitoggletool4,'State','off');
guidata(hObject,handles)

% Get main Database directory 
[PATH,Pt_name,~]=fileparts(fileparts(handles.fname));

% Create a JPEGs folder in the databse directory 
snap_dirname=char(fullfile(fileparts(PATH),'JPEGs'));
if ~exist(snap_dirname,'dir');  mkdir(snap_dirname); end

% Sets a default snapshot filename
snap_filename=char(fullfile(snap_dirname,strcat(Pt_name,'_Statistics.jpg')));

[filename,jpeg_dirname]=uiputfile({'*.jpg;','Image Files';...
          '*.*','All Files' },'Save Image',...
          snap_filename);

 snap_filename=char(fullfile(jpeg_dirname,filename));     

% If save operation is cancelled
if snap_filename(1) == 0 ;  return; end 
    
mia_export_fig(snap_filename,'-jpeg',handles.figure1)

guidata(hObject, handles);


% --------------------------------------------------------------------
function uitoggletool_displight_OnCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_displight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Reduce main figure in height and shift down
set(handles.figure1,'units','pixels');
% set(handles.axes_orig,'units', 'pixels');
% set(handles.axes_stats,'units', 'pixels');
% set(handles.axes_duration,'units', 'pixels');

fig_position = get(handles.figure1,'Position') ;
shift = fig_position(3) /3;
fig_position(3) = fig_position(3) - 2*shift ;
fig_position(1) = fig_position(1) + 2*shift ;
set(handles.figure1,'Position',fig_position) ;

pos = get(handles.axes_duration,'Position') ;
pos_orig = get(handles.axes_orig,'Position') ;

handles.backupAxesPos = get(handles.axes_duration,'Position');
% Set axes position (same x than first column in intial figure)
set(handles.axes_duration,'Position',[pos_orig(1)*3 pos(2) pos(3)*3 pos(4)]) ;

% Remove the two first axes and colorbars
set(handles.axes_orig,'visible','off') %hide the current axes
set(get(handles.axes_orig,'children'),'visible','off') ;%hide the current axes contents
set(handles.axes_stats,'visible','off') ;%hide the current axes
set(get(handles.axes_stats,'children'),'visible','off'); %hide the current axes contents
colorbar('peer',handles.axes_stats,'off');
colorbar('peer',handles.axes_orig,'off');

% Remove button to select chan (for now because it is not properly implemented)
set(handles.Select_chan, 'visible','off');

% Remove text (zscore + ttest uncorrected)
set(handles.text12, 'visible','off');
set(handles.text11, 'visible','off');
set(handles.edit_scale, 'visible','off');
set(handles.text_scalezs, 'visible','off');
set(handles.cursor_infos, 'visible','off');


% Put back the figure position mode to normalized 
set(handles.figure1,'units','normalized');

guidata(hObject, handles);


% --------------------------------------------------------------------
function uitoggletool_displight_OffCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_displight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Enlarge main figure in height and shift down
set(handles.figure1,'units','pixels');

fig_position = get(handles.figure1,'Position') ;
shift = fig_position(3) ;
fig_position(3) = fig_position(3) + 2*shift ;
fig_position(1) = fig_position(1) - 2*shift ;
set(handles.figure1,'Position',fig_position) ;

  % Set axes position (same x than first column in intial figure)
set(handles.axes_duration,'Position',handles.backupAxesPos) ;

% Make the two left axes and colorbars visible again
 set(handles.axes_orig,'visible','on') ;
 set(get(handles.axes_orig,'children'),'visible','on') ;
 set(handles.axes_stats,'visible','on') ;
 set(get(handles.axes_stats,'children'),'visible','on'); 
 colorbar('peer',handles.axes_orig,'location', 'NorthOutside');
 colorbar('peer',handles.axes_stats,'location', 'NorthOutside');
%  colorbar('peer',handles.axes_duration ,'location', 'NorthOutside');

%  set(handles.figure1,'units','normalized');
set(handles.figure1,'units','characters');

% Set Select button back to visibble 
 set(handles.Select_chan, 'visible','on');

% Set text back to visible (zscore + ttest uncorrected)
set(handles.text12, 'visible','on');
set(handles.text11, 'visible','on');
set(handles.edit_scale, 'visible','on');
set(handles.text_scalezs, 'visible','on');
set(handles.cursor_infos, 'visible','on');


% --- Executes on button press in pushbutton_display_timeseries.
function pushbutton_display_timeseries_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_display_timeseries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get all parameters for visulalization (from edit fields)

[threshp, threshdur,max_tvals,max_zscore] = get_parameters(handles) ; 

time = handles.Time; 
meanzs = handles.meanzsnew(handles.iSel,:,:);

set(0,'DefaultAxesColorOrder',jet(size(meanzs,1)));
 
% Plot timeseries of the selected channels
figure ; plot(time,meanzs','LineWidth',1.05); legend(strrep(handles.Labels(handles.iSel),'_','-')) ; grid on ; ylim([-max_zscore max_zscore]) ; 
xlabel('Time (sec)'); ylabel('Z-score');


function mouseMove (object, eventdata, handles)

% Get cursor position in the window
C = get (gcf, 'CurrentPoint');

% Set axes positions to pixels in order to get mouse cursor coordinates
set(handles.axes_orig,'units', 'characters');
set(handles.axes_stats,'units', 'characters');
set(handles.axes_duration,'units', 'characters');

% Get position of the axes 
xOrig = get(handles.axes_orig,'Position');
xStats = get(handles.axes_stats,'Position');
xDur = get(handles.axes_duration,'Position');

% Set axes positions back to normalized (auto adjust size with size window)
set(handles.axes_orig,'units', 'normalized');
set(handles.axes_stats,'units', 'normalized');
set(handles.axes_duration,'units', 'normalized');

% Get info on one of the panels (axes_orig)
xOrigLim = get(handles.axes_orig,'xlim');

% Duration in time
dur = sum(abs(xOrigLim)) ;

% Starting position of each panel
starts = [xOrig(1),xStats(1),xDur(1)];

% Cursor if out of the window
if find(starts-C(1)>0,1)-1 ==0 ; return; end

% Position of the cursor inside the panel
pos = C(1,1) - starts(find(starts-C(1)>0,1)-1);

if isempty(find(starts-C(1)>0,1)) 
    pos = C(1,1) - starts(3) ; 
end

% Lastency (in msec)
lat = pos*dur/xDur(3) + xOrigLim(1); 

% If mouse is not moving over one panel
if lat<xOrigLim(1)|lat>xOrigLim(2)
    set(handles.cursor_infos,'String',sprintf('Latency = ------'));
else
    % Display latency 
    set(handles.cursor_infos,'String',sprintf('Latency = %0.3f sec', lat));

end
