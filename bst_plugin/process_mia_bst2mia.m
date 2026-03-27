function varargout = process_mia_bst2mia(varargin)
% process_mia_bst2mia: Run mia_bst2mia from the Brainstorm pipeline.
% Brainstorm pipline which converts BST files to MIA format for the analysis using MIA.
%Inputs:
%   Condition        - String, condition name
%   ProtocolName     - String, Brainstorm protocol name 
%   LabelingTable    - String, full path to TSV labeling table
%   GroupChannelFile - String, full path to group channel MAT file
%   SubjectsToSkip   - Cell array of subject names to skip
%

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
% GetDescription: Define the Brainstorm process options and UI elements
% This function configures how the plugin appears and functions within Brainstorm's process pipeline
function sProcess = GetDescription() %#ok<DEFNU>
    global GlobalData;  % Access Brainstorm's global database for protocol information

    % === BASIC PROCESS METADATA ===
    sProcess.Comment     = 'MIA: Convert from BST to MIA';  % Display name in Brainstorm menu
    sProcess.Category    = 'Custom';                       % Category in process list
    sProcess.SubGroup    = 'Test';                          % Sub-category grouping
    sProcess.Index       = 910;                             % Menu index/ordering
    sProcess.Description = 'https://github.com/MIA-iEEG/mia';  % Link to documentation

    % === PROCESS I/O CONFIGURATION ===
    sProcess.InputTypes  = {'import'};      % Input type: import data from Brainstorm protocol
    sProcess.OutputTypes = {'matrix'};      % Output type: matrix data structure
    sProcess.nInputs     = 1;               % Number of inputs required
    sProcess.nMinFiles   = 0;               % Minimum files required (0 = process uses protocol data)

    % === RETRIEVE AVAILABLE PROTOCOLS FROM BRAINSTORM ===
    % Query the Brainstorm global database to get all available protocol names
    if isfield(GlobalData, 'DataBase') && ...
       isfield(GlobalData.DataBase, 'ProtocolInfo') && ...
       ~isempty(GlobalData.DataBase.ProtocolInfo)
        % Extract protocol names from the database
        ProtocolNames = {GlobalData.DataBase.ProtocolInfo.Comment};
    else
        % Fallback if no protocols are available
        ProtocolNames = {'No protocol found'};
    end

    % === PROTOCOL SELECTION DROPDOWN ===
    % Create a combobox to let users select which Brainstorm protocol to process
    sProcess.options.protocolsel.Comment = 'Protocol:';           % Label for the dropdown
    sProcess.options.protocolsel.Type    = 'combobox';            % UI element type
    sProcess.options.protocolsel.Value   = {1, ProtocolNames};    % Default: first protocol

    % === CONDITION NAME INPUT ===
    % Text field for specifying which experimental condition to process
    sProcess.options.condition.Comment = 'Condition name: ';   % Label
    sProcess.options.condition.Type    = 'text';              % UI element type
    sProcess.options.condition.Value   = '';                  % Default: empty (user must enter)

    % === TSV LABELING TABLE FILE BROWSER ===
    % File browser to select the channel/electrode labeling table
    SelectTSV = { ...
    '', ...                  % Default path (empty = current directory)
    'files', ...             % Selection mode: files only
    { ...
        '(*.tsv) TSV files (*.tsv)', '*.tsv'; ...    % File filter for TSV files
        'All files (*.*)', '*.*' ...                 % Allow any file type
    }, ...
    'DataIn'};               % Starting directory (data input folder)

    sProcess.options.tsvfile.Comment = 'Labeling table (TSV): ';  % Label
    sProcess.options.tsvfile.Type    = 'filename';               % UI element type
    sProcess.options.tsvfile.Value   = SelectTSV;                % File browser configuration

    % === GROUP CHANNEL FILE BROWSER ===
    % File browser to select the MAT file containing group electrode location data
    SelectChan = { ...
    '', ...                  % Default path (empty = current directory)
    'files', ...             % Selection mode: files only
    { ...
        '(*.mat) MAT files (*.mat)', '*.mat'; ...    % File filter for MAT files
        'All files (*.*)', '*.*' ...                 % Allow any file type
    }, ...
    'DataIn'};               % Starting directory (data input folder)

    sProcess.options.chanfile.Comment = 'Group channel file (MAT): ';  % Label
    sProcess.options.chanfile.Type    = 'filename';                   % UI element type
    sProcess.options.chanfile.Value   = SelectChan;                   % File browser configuration

    % === SUBJECT EXCLUSION LIST INPUT ===
    % Text field for specifying subjects to skip during conversion
    sProcess.options.excludesubs.Comment = ...
        'Subjects to exclude (comma-separated, leave blank = use all): ';  % Label with usage hint
    sProcess.options.excludesubs.Type  = 'text';     % UI element type
    sProcess.options.excludesubs.Value = '';         % Default: empty (process all subjects)
end


%% ===== FORMAT COMMENT =====
% FormatComment: Return the process display name
% This function is called by Brainstorm to get the text displayed in the process list
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    % Simply return the process comment/name defined in GetDescription
    Comment = sProcess.Comment;
end


%% ===== RUN =====
% Run: Main execution function called by Brainstorm pipeline
% Reads user inputs, validates them, and executes the MIA conversion
function OutputFiles = Run(sProcess, ~) %#ok<DEFNU>
    OutputFiles = {};  % Initialize empty output (will be populated if successful)

    % === STEP 1: EXTRACT USER-SELECTED OPTIONS FROM BRAINSTORM GUI ===
    % Parse protocol selection from the combobox (may be cell array)
    protocolOpt = sProcess.options.protocolsel.Value;
    if iscell(protocolOpt)
        % If it's a cell array, extract using index: {index, options}, get options[index]
        ProtocolName = strtrim(protocolOpt{2}{protocolOpt{1}});
    else
        % If it's a string, use directly
        ProtocolName = strtrim(protocolOpt);
    end

    % Extract other parameters and trim whitespace
    Condition  = strtrim(sProcess.options.condition.Value);        % Condition name
    TSVFile    = sProcess.options.tsvfile.Value{1};                % Path to TSV labeling table
    ChanFile   = sProcess.options.chanfile.Value{1};               % Path to group channel MAT file
    ExcludeRaw = strtrim(sProcess.options.excludesubs.Value);      % Comma-separated subject list

    % === STEP 2: VALIDATE ALL MANDATORY INPUTS ===
    % Check that required parameters are provided and accessible
    if isempty(ProtocolName)
        bst_report('Error', sProcess, [], 'Protocol name is empty.');
        return  % Exit early if validation fails
    end
    if isempty(Condition)
        bst_report('Error', sProcess, [], 'Condition name is empty.');
        return  % Exit early if validation fails
    end
    if isempty(TSVFile) || ~exist(TSVFile, 'file')
        bst_report('Error', sProcess, [], 'TSV labeling table not found.');
        return  % Exit early if validation fails
    end
    if isempty(ChanFile) || ~exist(ChanFile, 'file')
        bst_report('Error', sProcess, [], 'Group channel file not found.');
        return  % Exit early if validation fails
    end

    % === STEP 3: PARSE SUBJECT EXCLUSION LIST ===
    % Convert comma-separated string into a cell array for processing
    if isempty(ExcludeRaw)
        % If no subjects specified, process all subjects
        SubjectsToSkip = {};
    else
        % Split on commas and remove whitespace from each subject name
        SubjectsToSkip = strtrim(strsplit(ExcludeRaw, ','));
        % Remove any empty entries that may have been created
        SubjectsToSkip(cellfun(@isempty, SubjectsToSkip)) = [];
    end

    % Log the subjects being skipped to Brainstorm report
    bst_report('Info', sProcess, [], ...
        sprintf('Subjects to skip: %s', strjoin(SubjectsToSkip, ', ')));

    % === STEP 4: CALL MIA CONVERSION FUNCTION ===
    % Execute the main conversion with try-catch to handle errors gracefully
    % Function signature: mia_bst2mia(Condition, ProtocolName, LabelingTable, GroupChannelFile, SubjectsToSkip)
    try
        % Call the MIA conversion function and capture output structures
        [rois, groupChanTable, groupFdata, stat] = mia_bst2mia( ...
            Condition, ...           % Experimental condition to process
            ProtocolName, ...        % Brainstorm protocol name
            TSVFile, ...             % Path to labeling/electrode table
            ChanFile, ...            % Path to group channel file
            SubjectsToSkip);         % Subjects to exclude

        % === DISPLAY SUCCESS AND RESULTS ===
        bst_report('Info', sProcess, [], ...
            sprintf('Successfully processed condition: %s', Condition));
        disp('ROI structure created:');  % Display ROI data to MATLAB console
        disp(rois);

    catch ME
        % === HANDLE ERRORS ===
        % Catch any exceptions and report to Brainstorm
        bst_report('Error', sProcess, [], ...
            ['mia_bst2mia failed: ' ME.message]);
        return  % Exit without creating output files
    end
end