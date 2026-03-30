function varargout = process_mia_group_gui(varargin)
% process_mia_group_gui: Launch mia_group_gui from Brainstorm ROI files.

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'MIA: Visualize Averages';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Test';
    sProcess.Index       = 911;
    sProcess.Description = 'https://github.com/MIA-iEEG/mia';

    % No Process1 input is required: ROI files are discovered from the
    % current Brainstorm protocol and the selected subject.
    sProcess.InputTypes  = {'import'};
    sProcess.OutputTypes = {'import'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 0;

    [subjectNames, iDefaultSubject] = get_protocol_subject_options();
    sProcess.options.roisubject.Comment = 'ROI subject: ';
    sProcess.options.roisubject.Type    = 'combobox';
    sProcess.options.roisubject.Value   = {iDefaultSubject, subjectNames};

    sProcess.options.label1.Comment = [ ...
        '<BR>Select the subject whose ROI files should be visualized. <BR>' ...
        'After clicking Run, the process scans data/<subject>/ROIS and opens a checkbox dialog for the available *_rois.mat files. <BR>'];
    sProcess.options.label1.Type = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<INUSD,DEFNU>
    OutputFiles = {};

    % Read the currently active Brainstorm protocol to locate the subject
    % ROI folders on disk.
    prot = bst_get('ProtocolInfo');
    if isempty(prot) || ~isfield(prot, 'STUDIES') || isempty(prot.STUDIES)
        bst_report('Error', sProcess, [], 'No active Brainstorm protocol found.');
        return
    end

    % Extract the subject selected in the dropdown. Brainstorm stores
    % combobox values as {selectedIndex, optionList}.
    subjectOpt = sProcess.options.roisubject.Value;
    if iscell(subjectOpt) && numel(subjectOpt) >= 2 && ~isempty(subjectOpt{2})
        subjectNames = subjectOpt{2};
        if isempty(subjectOpt{1}) || subjectOpt{1} < 1 || subjectOpt{1} > numel(subjectNames)
            bst_report('Error', sProcess, [], 'Invalid ROI subject selection.');
            return
        end
        SubjectName = strtrim(subjectNames{subjectOpt{1}});
    else
        SubjectName = strtrim(subjectOpt);
    end

    if isempty(SubjectName)
        bst_report('Error', sProcess, [], 'ROI subject is empty.');
        return
    end

    % ROI files are stored in data/<subject>/ROIS.
    roiDir = fullfile(prot.STUDIES, SubjectName, 'ROIS');
    if ~exist(roiDir, 'dir')
        bst_report('Error', sProcess, [], ...
            sprintf('ROIS folder not found for subject "%s": %s', SubjectName, roiDir));
        return
    end

    [roiFiles, conditionNames] = get_subject_roi_files(roiDir);
    if isempty(roiFiles)
        bst_report('Error', sProcess, [], ...
            sprintf('No *_rois.mat files were found in %s', roiDir));
        return
    end

    % Let the user choose which ROI conditions to visualize.
    selectedIdx = select_roi_conditions(conditionNames, SubjectName);
    if isempty(selectedIdx)
        bst_report('Info', sProcess, [], 'No ROI conditions were selected.');
        return
    end

    selectedConditionNames = conditionNames(selectedIdx);
    selectedRois = cell(1, numel(selectedIdx));
    selectedPaths = cell(1, numel(selectedIdx));

    for iSel = 1:numel(selectedIdx)
        selectedPaths{iSel} = fullfile(roiDir, roiFiles(selectedIdx(iSel)).name);
        % Each ROI file is expected to contain a variable named "rois".
        roiData = load(selectedPaths{iSel}, 'rois');
        if ~isfield(roiData, 'rois')
            bst_report('Error', sProcess, [], ...
                sprintf('File does not contain a "rois" variable: %s', selectedPaths{iSel}));
            return
        end
        selectedRois{iSel} = roiData.rois;
    end

    assignin('base', 'mia_visualize_subject', SubjectName);
    assignin('base', 'mia_visualize_roi_paths', selectedPaths);
    assignin('base', 'mia_visualize_conditions', selectedConditionNames);

    % Match the mia_group_gui API:
    %   single condition  -> mia_group_gui(rois, 'Condition')
    %   multiple conditions -> mia_group_gui(rois1, rois2, ..., 'Cond1-Cond2')
    if numel(selectedRois) == 1
        mia_group_gui(selectedRois{1}, selectedConditionNames{1});
    else
        mia_group_gui(selectedRois{:}, strjoin(selectedConditionNames, '-'));
    end

    bst_report('Info', sProcess, [], ...
        sprintf('Opened MIA visualization for subject "%s": %s', ...
        SubjectName, strjoin(selectedConditionNames, ', ')));
end


% Return all protocol subjects so the user can pick which subject ROI
% folder should be scanned.
function [subjectNames, iDefaultSubject] = get_protocol_subject_options()
    subjectNames = {'No subject found'};
    iDefaultSubject = 1;

    try
        % Read available protocol subjects to populate the dropdown.
        subs = bst_get('ProtocolSubjects');
        if ~isempty(subs) && isfield(subs, 'Subject') && ~isempty(subs.Subject)
            subjectNames = {subs.Subject.Name};
            % Exclude Brainstorm group-analysis nodes from the subject list.
            subjectNames = subjectNames(~contains(subjectNames, 'Group'));
            if isempty(subjectNames)
                subjectNames = {'No subject found'};
            end
        end
    catch
        subjectNames = {'No subject found'};
    end
end


% Find all ROI files in data/<subject>/ROIS and expose the condition names
% by stripping the trailing "_rois.mat" suffix from each filename.
function [roiFiles, conditionNames] = get_subject_roi_files(roiDir)
    roiFiles = dir(fullfile(roiDir, '*_rois.mat'));
    roiFiles = roiFiles(~[roiFiles.isdir]);

    if isempty(roiFiles)
        conditionNames = {};
        return
    end

    [~, sortIdx] = sort(lower({roiFiles.name}));
    roiFiles = roiFiles(sortIdx);
    conditionNames = cell(size(roiFiles));
    for iFile = 1:numel(roiFiles)
        % Convert filenames like "Ap_bipolar_2_rois.mat" to "Ap_bipolar_2".
        conditionNames{iFile} = regexprep(roiFiles(iFile).name, '_rois\.mat$', '');
    end
end


% Open a small checkbox dialog so the user can choose one or more ROI
% conditions to visualize with mia_group_gui.
function selectedIdx = select_roi_conditions(conditionNames, subjectName)
    selectedIdx = [];

    nConditions = numel(conditionNames);
    if nConditions == 0
        return
    end

    dlgHeight = min(max(170 + 28 * nConditions, 260), 700);
    dlg = dialog( ...
        'Name', sprintf('Select ROI conditions: %s', subjectName), ...
        'Units', 'pixels', ...
        'Position', [200, 120, 420, dlgHeight], ...
        'WindowStyle', 'normal');

    uicontrol( ...
        'Parent', dlg, ...
        'Style', 'text', ...
        'String', 'Select the ROI conditions to visualize:', ...
        'HorizontalAlignment', 'left', ...
        'Units', 'normalized', ...
        'Position', [0.07, 0.90, 0.86, 0.06]);

    checkboxHandles = gobjects(nConditions, 1);
    topY = 0.84;
    rowStep = min(0.72 / max(nConditions, 1), 0.065);

    for iCond = 1:nConditions
        yPos = topY - (iCond - 1) * rowStep;
        checkboxHandles(iCond) = uicontrol( ...
            'Parent', dlg, ...
            'Style', 'checkbox', ...
            'String', conditionNames{iCond}, ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'left', ...
            'Position', [0.08, yPos, 0.84, 0.06], ...
            'Value', iCond == 1);
    end

    uicontrol( ...
        'Parent', dlg, ...
        'Style', 'pushbutton', ...
        'String', 'Cancel', ...
        'Units', 'normalized', ...
        'Position', [0.52, 0.05, 0.18, 0.08], ...
        'Callback', @(~, ~) cancel_selection());

    uicontrol( ...
        'Parent', dlg, ...
        'Style', 'pushbutton', ...
        'String', 'Visualize', ...
        'Units', 'normalized', ...
        'Position', [0.72, 0.05, 0.18, 0.08], ...
        'Callback', @(~, ~) confirm_selection());

    uiwait(dlg);

    if isvalid(dlg)
        if isappdata(dlg, 'selectedIdx')
            selectedIdx = getappdata(dlg, 'selectedIdx');
        end
        delete(dlg);
    end

    function confirm_selection()
        % Collect the indices of all checked conditions and return them to
        % the caller when the dialog closes.
        selectedIdx = find(arrayfun(@(h) logical(get(h, 'Value')), checkboxHandles));
        setappdata(dlg, 'selectedIdx', selectedIdx);
        uiresume(dlg);
    end

    function cancel_selection()
        % Return an empty selection when the user cancels.
        setappdata(dlg, 'selectedIdx', []);
        uiresume(dlg);
    end
end
