% LOAD_RIVM_DATA Loads the data from the RIVM-api's. 
%    It will load the data from the cache if it exists. If not, it 
%    will download the data using the `retrieve_rivm_data`-script.

disp('LOADING THE RIVM DATA...');
disp('  ');

%% Retrieve the data.
cacheFiles = getCacheFiles();
if ~isempty(cacheFiles)
    cacheFile = table2struct(cacheFiles(1,:));
    disp(['Loading data from file "' cacheFile.folder '/' cacheFile.name '".']);
    load([cacheFile.folder '/' cacheFile.name]);
else
    disp('No cache files found -> Downloading data from rivm...');
    disp('   ');
    scr.retrieve_rivm_data;
end

%% Remove clutter from the workspace
clear('cacheFiles');


%% Helper Functions

function result = getCacheFiles()
    if ~~exist('cache/rivm','dir')
        items = dir('cache/rivm');
        cacheFiles = items(arrayfun(@isCacheFile, items));
        if ~isempty(cacheFiles)
            result = sortrows(struct2table(cacheFiles), 'datenum', 'descend');
        else
            result = [];
        end
    else
        result = [];
    end
    
end


function result = isCacheFile(x)
    % Check if it is a mat file.
    if ~x.isdir && strcmp(x.name(end-2:end), 'mat')
        % Check if the file contains the nescesary variables.
        items = struct2table(whos('-file', [x.folder '/' x.name]));
        hasVars = [
            ismember('rivm_cases', items.name)
            ismember('rivm_infectiousPeople', items.name)
            ismember('rivm_reproduction', items.name)
            ismember('rivm_sewage', items.name)
            ismember('rivm_totalPerMunicipality', items.name)
        ];
        result = all(hasVars);
    else
        result = false;
    end
end