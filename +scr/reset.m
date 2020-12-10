% RESET Resets all the tasks that the `init` script set up.

%% Clear CBS Data
scr.clear_cbs_data;

%% Clear RIVM Data
scr.clear_rivm_data;

%% Clear Workspace
disp('Clearing the workspace...');
clear();