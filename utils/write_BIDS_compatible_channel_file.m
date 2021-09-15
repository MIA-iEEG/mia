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
% 2021/9/14 : ASD creation

function output = write_BIDS_compatible_channel_file(labels, isGood, fname, tag)

% Prepare list of labels and corresponding status (BIDS compatibility)) 
name = labels ; 
status_description = repmat({'n/a'},length(name),1);
status_description(~isGood) = {tag} ; 

% Create a two colum table (name, status)
T = table(name',status_description);

% Prepare output filename
[dir_name, pt_fname] = fileparts(fname);
output = fullfile(dir_name,strrep(pt_fname,'signal_LFP', 'channels.tsv')) ; 

% Write table <SUBJ>_channels.tsv
writetable(T,output, 'FileType','text','Delimiter','\t');
