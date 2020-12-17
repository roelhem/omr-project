% LOAD_RIVM_DATA Loads the data from the RIVM-api's. 
%    It will load the data from the cache if it exists. If not, it 
%    will download the data using the `retrieve_rivm_data`-script.
global cbs_populationTotal;
global cbs_AgeGroupPopulation;
disp('LOADING THE CBS DATA...');
disp('  ');

%% Retrieve the data.
cacheFiles = getCacheFiles();
if ~isempty(cacheFiles)
    cacheFile = table2struct(cacheFiles(1,:));
    disp(['Loading data from file "' cacheFile.folder '/' cacheFile.name '".']);
    load([cacheFile.folder '/' cacheFile.name]);
else
    disp('No cache files found -> Downloading data from cbs...');
    disp('   ');
    scr.retrieve_cbs_data;
end

%% Create the default classes.
cbs_AgeGroupPopulation = lib.classes.AgeGroupPopulation(...
    cbs_populationTotal.One.Group.Boundaries,... Loads the age group boundaries.
    cbs_populationTotal.One.Total...             Loads the total per age group.
);

%% Remove clutter from the workspace
clear('cacheFiles', 'cacheFile');


%% Helper Functions

function result = getCacheFiles()
    if ~~exist('cache/cbs','dir')
        items = dir('cache/cbs');
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
            ismember('cbs_populationTotal', items.name)
        ];
        result = all(hasVars);
    else
        result = false;
    end
end