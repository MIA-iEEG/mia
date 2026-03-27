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
% Subjects to exclude can be provided as a comma-separated string in the GUI (e.g. 'SUBJ1, SUBJ2').

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    global GlobalData;

    sProcess.Comment     = 'MIA: Convert from BST to MIA';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Test';
    sProcess.Index       = 910;
    sProcess.Description = 'https://github.com/MIA-iEEG/mia';

    sProcess.InputTypes  = {'import'};
    sProcess.OutputTypes = {'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 0;

    % === GET PROTOCOL LIST
    if isfield(GlobalData, 'DataBase') && ...
       isfield(GlobalData.DataBase, 'ProtocolInfo') && ...
       ~isempty(GlobalData.DataBase.ProtocolInfo)

        ProtocolNames = {GlobalData.DataBase.ProtocolInfo.Comment};
    else
        ProtocolNames = {'No protocol found'};
    end

    % === PROTOCOL SELECTION
    sProcess.options.protocolsel.Comment = 'Protocol:';
    sProcess.options.protocolsel.Type    = 'combobox';
    sProcess.options.protocolsel.Value   = {1, ProtocolNames};

    % === CONDITION NAME
    sProcess.options.condition.Comment = 'Condition name: ';
    sProcess.options.condition.Type    = 'text';
    sProcess.options.condition.Value   = '';

    % === TSV FILE
    SelectTSV = { ...
    '', ...
    'files', ...
    { ...
        '(*.tsv) TSV files (*.tsv)', '*.tsv'; ...  % Note: semicolon separates name/pattern
        'All files (*.*)', '*.*' ...
    }, ...
    'DataIn'};

    sProcess.options.tsvfile.Comment = 'Labeling table (TSV): ';
    sProcess.options.tsvfile.Type    = 'filename';
    sProcess.options.tsvfile.Value   = SelectTSV;

    % === CHANNEL FILE
    SelectChan = { ...
    '', ...
    'files', ...
    { ...
        '(*.mat) MAT files (*.mat)', '*.mat'; ...  % Note: semicolon separates name/pattern
        'All files (*.*)', '*.*' ...
    }, ...
    'DataIn'};

    sProcess.options.chanfile.Comment = 'Group channel file (MAT): ';
    sProcess.options.chanfile.Type    = 'filename';
    sProcess.options.chanfile.Value   = SelectChan;

    % === SUBJECT EXCLUSION LIST
    sProcess.options.excludesubs.Comment = ...
        'Subjects to exclude (comma-separated, leave blank = use all): ';
    sProcess.options.excludesubs.Type  = 'text';
    sProcess.options.excludesubs.Value = '';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, ~) %#ok<DEFNU>
    OutputFiles = {};

    % 1. Read options
    protocolOpt = sProcess.options.protocolsel.Value;
    if iscell(protocolOpt)
        ProtocolName = strtrim(protocolOpt{2}{protocolOpt{1}});
    else
        ProtocolName = strtrim(protocolOpt);
    end

    Condition  = strtrim(sProcess.options.condition.Value);
    TSVFile    = sProcess.options.tsvfile.Value{1};
    ChanFile   = sProcess.options.chanfile.Value{1};
    ExcludeRaw = strtrim(sProcess.options.excludesubs.Value);

    % 2. Validate mandatory inputs
    if isempty(ProtocolName)
        bst_report('Error', sProcess, [], 'Protocol name is empty.');
        return
    end
    if isempty(Condition)
        bst_report('Error', sProcess, [], 'Condition name is empty.');
        return
    end
    if isempty(TSVFile) || ~exist(TSVFile, 'file')
        bst_report('Error', sProcess, [], 'TSV labeling table not found.');
        return
    end
    if isempty(ChanFile) || ~exist(ChanFile, 'file')
        bst_report('Error', sProcess, [], 'Group channel file not found.');
        return
    end

    % 3. Convert reject string to cell array
    if isempty(ExcludeRaw)
        SubjectsToSkip = {};
    else
        SubjectsToSkip = strtrim(strsplit(ExcludeRaw, ','));
        SubjectsToSkip(cellfun(@isempty, SubjectsToSkip)) = [];
    end

    bst_report('Info', sProcess, [], ...
        sprintf('Subjects to skip: %s', strjoin(SubjectsToSkip, ', ')));

    % 4. Call mia_bst2mia with correct argument order
    % mia_bst2mia(Condition, ProtocolName, LabelingTable, GroupChannelFile, SubjectsToSkip)
    try
        [rois, groupChanTable, groupFdata, stat] = mia_bst2mia( ...
            Condition, ...
            ProtocolName, ...
            TSVFile, ...
            ChanFile, ...
            SubjectsToSkip);

        % Display results
        bst_report('Info', sProcess, [], ...
            sprintf('Successfully processed condition: %s', Condition));
        disp('ROI structure created:');
        disp(rois);

    catch ME
        bst_report('Error', sProcess, [], ...
            ['mia_bst2mia failed: ' ME.message]);
        return
    end
end