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
%
% Code inspired by nst_get_version (NIRSTORM brainstorm plugin) 

function version = mia_get_version()

%MIA_GET_VERSION Return the current version of MIA
version = 'github-master';

% Get MIA root dir
path = fileparts(which('mia'));

% Open VERSION file (root dir) 
id = fopen(fullfile(fileparts(path), 'VERSION'));

% Get version number (after stabel_) 
str=fread(id,'*char' )';
tmp = strsplit(str,'_') ; 
version = tmp{2};

end

