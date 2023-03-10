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
% Copyright (C) 2016-2022 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)

function [] = mia_fix_database(group_dir)

% Reads all folders that are in indir 
d = dir(group_dir); 
istudies = [d(:).isdir]; % returns logical vector if is folder
studies = {d(istudies).name}';
studies(ismember(studies,{'.','..'})) = []; % Removes . and ..

%Loop through studies
for ss=1:length(studies) 

    d = dir(fullfile(group_dir,studies{ss},strcat(studies{ss},'*.mat'))); 
    dat  = load(fullfile(group_dir,studies{ss},d.name)) ; 
    
    for gg=1:length(dat.ganalysis)
        
        % Get the name of the folder containing subjects 
        [~,db_name,~] = fileparts(fileparts(fileparts(dat.ganalysis{gg}.fname)));
        [~,wrong_subject,~] = fileparts(fileparts(dat.ganalysis{gg}.fname)) ; 
        correct_subject = dat.ganalysis{gg}.subj ; 
        
        [~,fname,~] = fileparts(dat.ganalysis{gg}.fname);
      
        dat.ganalysis{gg}.fname = fullfile(fileparts(group_dir),db_name,correct_subject, strrep(fname,wrong_subject,correct_subject));
      
    end
    
    % Save back study
    save(fullfile(group_dir,studies{ss},d.name),'-struct','dat');
    
end

