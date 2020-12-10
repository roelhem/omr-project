% CLEAR_CBS_DATA Clears all the caches containing the cbs data.

%% Remove the cache folder
if exist('cache/cbs', 'dir')
    disp('Removing the CBS-data caches...');
    rmdir('cache/cbs', 's');
else
    disp('CBS-data caches were already removed...');
end

%% Remove variables from the workspace
disp('Removing CBS data from the workspace...');
clear('cbs_populationTotal');