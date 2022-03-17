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
% mia verbose mode (inspired by eegh in EEGLAB)

function cmd_history = mia_cmd_history(cmd)

global MIA_HISTORY_CMD;

% Prints the history onto the Command Window
if nargin < 1
	if isempty(MIA_HISTORY_CMD)
		fprintf('No history\n');
	else	
      for index = 1:length(MIA_HISTORY_CMD)
         txt = MIA_HISTORY_CMD{ length(MIA_HISTORY_CMD)-index+1 };
      end
    end	
    cmd_history = char(MIA_HISTORY_CMD);
    
% Add something to the history
elseif nargin == 1
	if isempty(cmd)
		return;
	end
	if ischar(cmd)
    	if isempty(MIA_HISTORY_CMD)
			MIA_HISTORY_CMD = {sprintf('%s : %s', datetime,cmd)};
		else	
			MIA_HISTORY_CMD = {sprintf('%s : %s', datetime,cmd) MIA_HISTORY_CMD{:}};
        end	
    end
end
