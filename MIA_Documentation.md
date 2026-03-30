# Using MIA concatenate channel function

Path: `MIA/bst_plugin/process_concatenate_channels.m`

This process creates one new grand subject that contains the channels of the selected subjects in the current Brainstorm protocol. This makes the later ROI-based MIA analysis easier to run on a common channel space.

Brainstorm menu:
`Run -> Add process icon -> Standardize -> MIA: Concatenate Channels`

It contains 2 input fields:

1. **New subject name**
The custom name of the grand subject that will be created. Default name: `COREG`.

2. **Subjects to skip**
Comma-separated list of subjects that should not be included in the concatenation. Leave empty to include all subjects.

Example:
`Subject01, Subject02`

Then click `Run`.


# Using MIA: Convert from BST to MIA function

Path: `MIA/bst_plugin/process_mia_bst2mia.m`

This process converts Brainstorm data to MIA ROI data using the current Brainstorm protocol, the condition dropped in `Process1`, the selected subject channel file, and the labeling table.

Brainstorm menu:
`Run -> Add process icon -> Test -> MIA: Convert from BST to MIA`

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


