% CLEAR_RIVM_DATA Clears all the caches containing the rivm data.

%% Remove the cache folder
if exist('cache/rivm', 'dir')
    disp('Removing the RIVM-data caches...');
    rmdir('cache/rivm', 's');
else
    disp('RIVM-data caches were already removed...');
end

%% Remove variables from the workspace
disp('Removing RIVM data from the workspace...');
clear('rivm_cases',...
      'rivm_infectiousPeople',...
      'rivm_reproduction',...
      'rivm_sewage',...
      'rivm_totalPerMunicipality'...
);