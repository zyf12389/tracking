% This script can be used to pack the results and submit them to a challenge.

addpath('C:\Users\ww\Desktop\visual tracking\VOT\vot-toolkit-master\vot-toolkit-master'); toolkit_path; % Make sure that VOT toolkit is in the path

[sequences, experiments] = workspace_load();

tracker = tracker_load('NCC');

workspace_submit(tracker, sequences, experiments);

