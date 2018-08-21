% Import new data in MIA
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
% Copyright (C) 2016-2018 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)

function [handles] = mia_import_data(handles,~) 

% Open a browser window using BST functions 
inputFormat{1,1} = {'vhdr'};
inputFormat{1,2} = 'SEEG: BrainVision (*.vhdr)';

inputFormat{2,1} = {'*'};
inputFormat{2,2} = 'SEEG: Brainstorm directory (*.*)';

if isfield(handles,'history')&&isfield(handles.history,'datapath')
    datapath = handles.history.datapath;
else
    datapath = handles.extOPTIONS.outdir ; 
end

% BST function to open a Java file dialog box
[RawFiles, FileFormat] = dialog_import_data('MIA : Open SEEG recordings...', ...  % Window title
     datapath, ...           % Working directory
    inputFormat);    % List of available file formats
 
if isempty(RawFiles) ; return ; end 

% Brainvision single patient
if strcmp(FileFormat,inputFormat{1,2})

    for ff=1:size(RawFiles,1)

        % Get path and default subject (=filename)
        [PathName,defaultName,EXT] = fileparts(RawFiles{ff}) ;
        handles.history.datapath = PathName;

        % Check if .eeg exists
        if ~exist(fullfile(PathName,strcat(defaultName,'.eeg')),'file')
            errordlg('.eeg file is missing.');
            continue;
        end

        % Get the list of existing subjects
        subjects = get(handles.list_patient,'string');

        if ~isempty(subjects)
            % Ask for the patient name and check if it already exist or not
            Name=check_name(subjects,get(handles.outdir,'string'),defaultName);
        else 
            Name = defaultName; 
        end
        
        % Use of cancel button
        if isempty(Name) ; return; end

        % Create progress bar for patients processing
        hwait = waitbar(50,'','Units','Normalized','Name','Loading data....');

        % Extract data
        handles.extOPTIONS.maindir = PathName;
        handles.extOPTIONS.overwrite = 'Yes';
        [OutputFile] = mia_s1_extract_vhdr_data(RawFiles(ff),handles.extOPTIONS,Name);

        % Close progress bar
        delete(hwait) ;
    end
    
% Brainstorm single/multiple patients
elseif strcmp(FileFormat,inputFormat{2,2})

    for ff=1:length(RawFiles) 
    
        % Check the validity of the directory
        [status] = check_bst_dir(RawFiles{ff}); 

        if status ==1
            % Get the list of existing subjects
            subjects = get(handles.list_patient,'string');

            % Get the name of the directory where user browse
            [maindir,defaultName,~] = fileparts(RawFiles{ff}) ; 

            % Ask for the patient name and check if it already exist or not
            PtName=check_name(subjects,get(handles.outdir,'string'),defaultName);

            % Use of cancel button
            if isempty(PtName); return; end

            % Create progress bar for patients processing
            hwait = waitbar(50,'','Units','Normalized','Name','Loading data....');

            OPTIONS.maindir=RawFiles{ff};
            OPTIONS.outdir=get(handles.outdir,'String');
            OPTIONS.SensorType='SEEG';

            mia_s1_extract_bst_data(OPTIONS,PtName);

            % Close progress bar
            delete(hwait) ;
        end
    
    end
    
end

% --- Check the content of the brainstorm directory to import
function [status]=check_bst_dir(directoryname)
  
% list file in directory
d=dir(directoryname);
isub=[d(:).isdir];
files={d(~isub).name};
status = 1 ; 
 
if isempty(files)
    status =0 ; 
    errordlg('Invalid Folder');
elseif strcmpi(files,'channel.mat')==zeros(1,length(files))
    % Check if the channel.mat file exist
    status =-1 ; 
    errordlg('Channel file is missing');
elseif ~cellfun(@isempty,strfind(files,'trial'))==zeros(1,length(files))
    % Check if the trials exist
    status =-1 ; 
    errordlg('Trial files is missing');
end

% --- Check if patient ID exists
function [NAME]=check_name(NAMES,outdir,GivenName)

%get back the given name
NAME=inputdlg('Enter patient name:','PATIENT NAME',1,cellstr(GivenName));

%Use of cancel button
if isempty(NAME); return; end

%files existing in outdir
if ~isempty(NAMES)
    %if the name already exist :
    if sum(ismember(NAMES,NAME))==1
        
        % Creation of a list containing all the names that have been used in the
        %past
        list='';
        for ii=1:length(NAMES)
            list=[list, '  ',char(NAMES(ii))];
        end
        
        % ask for overwriting or choose an other one
        choice = questdlg(sprintf('Those names alerady exist :%s',list) , ...
            'PATIENT NAME','Overwrite','Choose an other Name','Overwrite');
        switch choice
            case 'Overwrite'
                %find back the full existing file name and delete it
                rmdir(char(fullfile(outdir,NAME)),'s');
                
            case 'Choose an other Name'
                %function to check if the given name already exist
                NAME=check_name(NAMES,outdir,GivenName);
        end
        
    end
end
