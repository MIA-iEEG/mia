# Using MIA concatenate channel function

Path: `MIA/bst_plugin/process_concatenate_channels.m`

This process creates one new grand subject that contains the channels of the selected subjects in the current Brainstorm protocol. This makes the later ROI-based MIA analysis easier to run on a common channel space.

Brainstorm menu:
`Run -> Add process icon -> Standardize -> MIA: Concatenate Channels`

<p float="left">
  <img src="Media/2 Grand subject.png" height="300" />
  <img src="Media/2_1.png" height="300" />
</p>

It contains 2 input fields:

1. **New subject name**
The custom name of the grand subject that will be created. Default name: `COREG`.

2. **Subjects to skip**
Comma-separated list of subjects that should not be included in the concatenation. Leave empty to include all subjects.

Example:
`Subject01, Subject02`

Then click `Run`.

# Calculating time/frequency with Morlet method

Path: `Mia/bst_plugin/process_mia_extract_tf.m`

This process calculates time/frequency representations of Brainstorm data using the Morlet wavelet method. It uses the current Brainstorm protocol, the condition dropped in `Process1`.

Brainstorm menu: `Frequency -> MIA: Time-frequency (Morlet by band + 1/f norm)`

<p float="left">
  <img src="Media/3 Calculating time frequency.png" height="300" />
  <img src="Media/3_1 .png" height="300" />
</p>

It contains the following input fields:

1. **Baseline**
The time window (in milliseconds) used as the baseline for 1/f normalization. You can specify a start and end time (e.g., `-400.0` to `-1.0` ms) or check the `All file` box to use the entire file as baseline.

2. **Frequency bands**
The frequency range and steps to extract. Follows the MATLAB array format `start:step:end` (e.g., `50:10:170` to extract frequencies from 50Hz to 170Hz with a 10Hz step).

3. **Number of cycles**
The number of cycles for the central frequency of the Morlet wavelet (e.g., `7`), which defines the time-frequency resolution trade-off.

Then click `Run`.


# Using MIA: Convert from BST to MIA function

Path: `MIA/bst_plugin/process_mia_bst2mia.m`

This process converts Brainstorm data to MIA ROI data using the current Brainstorm protocol, the condition dropped in `Process1`, the selected subject channel file, and the labeling table.

Brainstorm menu:
`Run -> Add process icon -> Test -> MIA: Convert from BST to MIA`

<p float="left">
  <img src="Media/4 Converting to MIA.png" height="300" />
  <img src="Media/4_1.png" height="300" />
</p>

## Current inputs

The process now contains 3 user-facing inputs:

1. **Files in Process1**
Drop one condition from the Brainstorm database. The process reads the condition automatically from the first dropped input with:
`sInputs(1).Condition`

1. **Labeling table (TSV)**
Path to the `.tsv` file containing the labeling information.

1. **Channel subject**
Dropdown listing the subjects in the current Brainstorm protocol. The selected subject is used to resolve the `channel.mat` file automatically from the protocol studies directory.

1. **Subjects to skip**
Comma-separated list of subjects to exclude from the conversion. Leave empty to use all eligible subjects.

Example:
`Subject01, Subject02`

Then click `Run`.


## Updated workflow in `process_mia_bst2mia.m`

The process now works as follows:

1. Read the current Brainstorm protocol automatically:
`prot = bst_get('ProtocolInfo');`

2. Read the condition automatically from the first dropped input:
`Condition = strtrim(sInputs(1).Condition);`

3. Read the selected channel subject from the dropdown.

4. Resolve the corresponding `channel.mat` automatically inside the current protocol studies directory.

5. Pass the resolved channel-file path to `mia_bst2mia(...)`.
   
6. The mia_bst2mia function then reads all subjects in the current protocol, excluding any specified in the "Subjects to skip" input, and converts their Brainstorm data for the selected condition into MIA ROI data and saves the resulting ROI fies in the 'data/group_subject(COREG)/ROIS' folder of the current protocol.
   
# Using MIA: Visualize Averages function

Path: `MIA/bst_plugin/process_mia_group_gui.m`

This process opens `mia_group_gui` directly from ROI files already saved in the current Brainstorm protocol.

Brainstorm menu:
`Run -> Add process icon -> Test -> MIA: Visualize Averages`

<p float="left">
  <img src="Media/5 Visualizing avg.png" height="200" />
  <img src="Media/5_1.png" height="200" />
  <img src="Media/5_2.png" height="200" />
  <img src="Media/5_3.png" height="200" />
  <img src="Media/5_4.png" height="200" />
  <img src="Media/5_5.png" height="200" />
</p>

## Current inputs

The process contains 1 user-facing input:

1. **ROI subject**
Dropdown listing the subjects in the current Brainstorm protocol. The selected subject is used to find ROI files inside:
`data/<SelectedSubject>/ROIS`

Then click `Run`.


## Updated workflow in `process_mia_group_gui.m`

The process now works as follows:

1. Read the current Brainstorm protocol automatically:
`prot = bst_get('ProtocolInfo');`

2. Read the selected subject from the dropdown.

3. Scan the ROI folder of the selected subject:
`data/<SelectedSubject>/ROIS`

4. Find all files matching:
`*_rois.mat`

5. Convert each filename into a condition name by removing the suffix:
`_rois.mat`

Example:
`Ap_bipolar_2_rois.mat -> Ap_bipolar_2`

6. Open a checkbox dialog so the user can select one or more conditions to visualize.

7. Load the `rois` variable from each selected file.

8. Call `mia_group_gui(...)` with the selected ROI structures.


## Behavior of the final call

If one condition is selected, the process calls:

```matlab
mia_group_gui(selectedRois{1}, selectedConditionNames{1});
```

Equivalent example:

```matlab
Ap_bipolar_2_file = load('data/COREG/ROIS/Ap_bipolar_2_rois.mat', 'rois');
mia_group_gui(Ap_bipolar_2_file.rois, 'Ap_bipolar_2');
```

If multiple conditions are selected, the process calls:

```matlab
mia_group_gui(selectedRois{:}, strjoin(selectedConditionNames, '-'));
```

Equivalent example:

```matlab
cond_1_file = load('data/COREG/ROIS/cond_1_rois.mat', 'rois');
cond_2_file = load('data/COREG/ROIS/cond_2_rois.mat', 'rois');
mia_group_gui(cond_1_file.rois, cond_2_file.rois, 'cond_1-cond_2');
```


## Helper functions added in `process_mia_group_gui.m`

### `get_protocol_subject_options()`

This helper reads all subjects from the current Brainstorm protocol and uses them to populate the `ROI subject` dropdown.


### `get_subject_roi_files(roiDir)`

This helper scans the selected subject ROI folder, finds all `*_rois.mat` files, sorts them, and extracts the condition names from the filenames.


### `select_roi_conditions(conditionNames, subjectName)`

This helper opens the checkbox dialog displayed after clicking `Run`. It lets the user choose which available ROI conditions should be passed to `mia_group_gui(...)`.


## Stats

need to install the packages:
1) signal processing toolbox (for filtfilt)
2) statistics and machine learning toolbox (for tinv)