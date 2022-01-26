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
% This function pops a dialog box to select files 
function [RawFiles, FileFormat] = mia_dialog_getfile(WindowTitle, DefaultDir, Filters)

import javax.swing.JFileChooser;

% Create the file selection figure 
hFig = figure('Name',WindowTitle, 'NumberTitle','off');

% Create a java file chooser object 
jfs = javaObjectEDT('com.mathworks.hg.util.dFileChooser');  % or: com.mathworks.hg.util.dFileChooser

% Displays it on the figure 
[hjFileChooser,hContainer] = javacomponent(jfs, [0,0,hFig.Position(3:4)], hFig);

% Add the approval button name
hjFileChooser.setApproveButtonText('Select'); 

% Set defautl directory 
hjFileChooser.setCurrentDirectory(java.io.File(DefaultDir));

% Allow multiple selection
hjFileChooser.setMultiSelectionEnabled(true);

% Removes 'All Files' from filter
hjFileChooser.setAcceptAllFileFilterUsed(false);

hjFileChooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);

for ff=1:length(Filters)
    
    fileFilter{ff} = com.mathworks.mwswing.FileExtensionFilter(Filters{ff,2}, Filters{ff,1}, false, true);
    hjFileChooser.addFileFilter(fileFilter{ff}) ; 

end

hjFileChooser.setFileFilter(fileFilter{1}) ;

hjFileChooser.ActionPerformedCallback = {@ActionPerformedCallback,hFig};

% % MAGIC: Print something to the console output to get it to flush something, if not it crashes on MacOS 10.9.2
% drawnow;
% fprintf(1, ' \b');

% Wait for user input
uiwait(hFig);
  
% We get here if the select button were pressed
if ishghandle(hFig)
  % Open were clicked
  RawFiles = getappdata(hFig,'selectedFile');
  FileFormat  = getappdata(hFig,'selectedFormat');
  close(hFig);
else  % figure was deleted/closed/canceled
  RawFiles = [] ;
  FileFormat  = [] ; 
end

end

%% Figure actions (Cancel & Open)
function ActionPerformedCallback(hjFileChooser, eventData, hFig)
  switch char(eventData.getActionCommand)
      case 'CancelSelection'
          close(hFig);

      case 'ApproveSelection'
          files = cellfun(@char, cell(hjFileChooser.getSelectedFiles), 'uniform',0);
          format = hjFileChooser.getFileFilter.getDescription ; 
          
          %msgbox(['Selected file: ' files{:}], 'Selected!');
%           if numel(files)==1
%               %files = files{1};
%           end
          if isempty(files)
              files = char(hjFileChooser.getSelectedFile);
          end
          setappdata(hFig,'selectedFile', files);
          setappdata(hFig,'selectedFormat', format);
          uiresume(hFig);

      otherwise
          % should never happen
  end
end  % ActionPerformedCallback

