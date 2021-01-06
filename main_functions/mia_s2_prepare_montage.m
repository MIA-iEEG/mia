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
function [OutputFile] = mia_s2_prepare_montage(varargin)

% Check if sInputs is empty and read the files for each patient
if nargin<2
    OPTIONS = varargin{1};
    
    % If subjects are specificed
    if isfield(OPTIONS,'subjects')
        subjects = OPTIONS.subjects;
    else % ALL subjects in directory
        d = dir(OPTIONS.outdir);
        isub = [d(:).isdir]; % returns logical vector if is folder
        subjects = {d(isub).name}';
        subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..
    end
    
    % Prepare input filename
    for ss=1:length(subjects) ;
        tmp = dir(fullfile(OPTIONS.outdir,subjects{ss},'*_signal_LFP.mat')) ;
        sInputs{ss} = fullfile(OPTIONS.outdir,subjects{ss},tmp.name) ;
    end
else % old but we keep for compatibility
    sInputs = varargin{1} ;
    OPTIONS = varargin{2};
end

% % Check if sInputs is empty and read the files for each patient
% if isempty (sInputs)
%     d = dir(OPTIONS.outdir); 
%     isub = [d(:).isdir]; % returns logical vector if is folder
%     subjects = {d(isub).name}';
%     subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..
%     for ss=1:length(subjects) ; 
%         file = strcat(subjects{ss},'_signal_LFP.mat');
%         tmp = fullfile(OPTIONS.outdir,subjects{ss},file) ; 
%         sInputs{ss} = tmp ; 
%     end
% end


% Loop through all subjects
for ii=1:length(sInputs)
    
    % Choose between Bipolar or Monopolar Montage    
    [PATHSTR,NAME,EXT] = fileparts(char(sInputs{ii})) ;
    
    % Prepare output filename : replace _signal_LFP by montage_LFP_data 
    outname = fullfile(PATHSTR,strrep(NAME,'_signal_LFP', strcat('_',OPTIONS.mtg,'_LFP_data'))) ;
  
    % Check if the files exist 
    if ((exist (strcat(outname, '.mat'), 'file')) && (strcmpi(OPTIONS.overwrite , 'Yes'))) ||(~exist (strcat(outname, '.mat'), 'file')) 
        % Load signal LFP
        load(sInputs{ii});
        fprintf('\nPreparing %s montage. Patient = %s  %d / %d\n\n',OPTIONS.mtg,NAME , ii, length(sInputs) ); 

        % Compute bipolar montage
        [dc, bilabels] = mia_make_bipolarmtg(F,labels);
        
        if strcmpi(OPTIONS.mtg,'Bipolar'); F = dc ; labels = bilabels ; end;

        % Save data 
        save(outname, 'Time', 'F','labels');
        fprintf('\nPatient = %s  %d / %d : %s montage prepared correctly\n\n',NAME, ii ,length(sInputs), OPTIONS.mtg );

    end
     OutputFile{ii} = outname ;   
end 
end
