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
    % === BASIC PROCESS METADATA ===
    sProcess.Comment     = 'MIA: Convert from BST to MIA';  % Display name in Brainstorm menu
    sProcess.Category    = 'Custom';                       % Category in process list
    sProcess.SubGroup    = 'Test';                          % Sub-category grouping
    sProcess.Index       = 910;                             % Menu index/ordering
    sProcess.Description = 'https://github.com/MIA-iEEG/mia';  % Link to documentation

    % === PROCESS I/O CONFIGURATION ===
    sProcess.InputTypes  = {'data'};        % Input type: recordings dropped in Process1
    sProcess.OutputTypes = {'matrix'};      % Output type: matrix data structure
    sProcess.nInputs     = 1;               % Number of inputs required
    sProcess.nMinFiles   = 1;               % Require at least one dropped file/folder in Process1

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

    % === SUBJECT USED TO RESOLVE GROUP CHANNEL FILE ===
    [subjectNames, iDefaultSubject] = get_protocol_subject_options();
    sProcess.options.channelsubject.Comment = 'Channel subject: ';
    sProcess.options.channelsubject.Type    = 'combobox';
    sProcess.options.channelsubject.Value   = {iDefaultSubject, subjectNames};

    % === SUBJECT EXCLUSION LIST INPUT ===
    % Text field for specifying subjects to skip during conversion
    sProcess.options.excludesubs.Comment = ...
        'Subjects to exclude (comma-separated, leave blank = use all): ';  % Label with usage hint
    sProcess.options.excludesubs.Type  = 'text';     % UI element type
    sProcess.options.excludesubs.Value = '';         % Default: empty (process all subjects)

    % === HELP TEXT AT THE BOTTOM
    sProcess.options.label1.Comment = [ ...
        '<BR>This process converts Brainstorm data to MIA format for analysis. <BR>' ...
        'The current Brainstorm protocol and the first dropped condition are used automatically. <BR>' ...
        'Select the subject whose channel.mat should be used, then provide the corresponding .tsv labelling table. <BR>' ...
        'Subjects to exclude can be listed separated by commas. <BR>' ...
    ];

    sProcess.options.label1.Type = 'label';

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
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};  % Initialize empty output (will be populated if successful)

    % === STEP 1: EXTRACT USER-SELECTED OPTIONS FROM BRAINSTORM GUI ===
    % Use the current Brainstorm protocol and the first dropped condition.
    prot = bst_get('ProtocolInfo');
    ProtocolName = strtrim(prot.Comment);
    Condition  = strtrim(sInputs(1).Condition);                    % Condition from dropped input
    subjectOpt = sProcess.options.channelsubject.Value;
    if iscell(subjectOpt)
        ChannelSubject = strtrim(subjectOpt{2}{subjectOpt{1}});
    else
        ChannelSubject = strtrim(subjectOpt);
    end
    TSVFile    = sProcess.options.tsvfile.Value{1};                % Path to TSV labeling table
    ChanFile   = resolve_subject_channel_file(prot, ChannelSubject);
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
    if isempty(ChannelSubject)
        bst_report('Error', sProcess, [], 'Channel subject is empty.');
        return
    end
    if isempty(TSVFile) || ~exist(TSVFile, 'file')
        bst_report('Error', sProcess, [], 'TSV labeling table not found.');
        return  % Exit early if validation fails
    end
    if isempty(ChanFile) || ~exist(ChanFile, 'file')
        bst_report('Error', sProcess, [], ...
            sprintf('Channel file not found for subject "%s".', ChannelSubject));
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


% Get all subjects from the current Brainstorm protocol to populate the
% dropdown used for channel-file selection. If subjects are found, default
% to the first one in the list.
function [subjectNames, iDefaultSubject] = get_protocol_subject_options()
    subjectNames = {};
    iDefaultSubject = [];

    try
        subs = bst_get('ProtocolSubjects');

        if ~isempty(subs) && isfield(subs, 'Subject') && ~isempty(subs.Subject)
            subjectNames = {subs.Subject.Name};
        end

        if ~isempty(subjectNames)
            iDefaultSubject = 1;
        end
    catch
        % Return empty values if no protocol is loaded or Brainstorm is not
        % ready to provide the subject list yet.
        subjectNames = {};
        iDefaultSubject = [];
    end
end

% Find the channel.mat file for the selected subject inside the current
% protocol studies directory. Search the most common Brainstorm study
% folders first and return the first match that exists.
function ChanFile = resolve_subject_channel_file(prot, subjectName)
    ChanFile = '';

    if isempty(subjectName) || isempty(prot) || ~isfield(prot, 'STUDIES') || isempty(prot.STUDIES)
        return
    end

    subjectStudyDir = fullfile(prot.STUDIES, subjectName);
    if ~exist(subjectStudyDir, 'dir')
        return
    end

    % Search order: raw matrix studies
    % then any direct child study folder containing channel.mat.
    searchPatterns = {fullfile(subjectStudyDir, '@rawmatrix*', 'channel.mat')};

    for iPattern = 1:numel(searchPatterns)
        channelMatches = dir(searchPatterns{iPattern});
        channelMatches = channelMatches(~[channelMatches.isdir]);
        if ~isempty(channelMatches)
            ChanFile = fullfile(channelMatches(1).folder, channelMatches(1).name);
            return
        end
    end
end
