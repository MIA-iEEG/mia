

# Using MIA concatenate channel function 

path: MIA/bst_plugin/process_concatenate_channels.m

- This functionis used to create one new subject which contains all the channels of the original subjects. This is done in oder for the MIA pipeline to adapt to analsyis ROI wise.

Run → Add process icon → Standardized → MIA: Concatenate Channels

It contains 2 input fields:

1) **New Subject name**: Which is the custom name that could be given to the new grand subject that will be created. Default name ```COREG```.
2) **Subjects to skip**: Here the user could specify the subjects that they want to exclude from the concatenation. The default is empty, meaning that all subjects will be included in the concatenation.

   - The subjects to skip can be written eg. as ```Subject01, Subject02```.

Then hit ```Run``` to execute.