
%% Loading the data
disp('Loading the data from the CBS API...');

disp('Loading population totals of the Netherlands...');
cbs_populationTotal = lib.loaders.cbs_populationTotal;

%% Determine the filename
dateName = datestr(datetime,"yyyy-mm-ddTHH-MM-ss");
dirName = 'cache/cbs/';
fileName = [dirName dateName '.mat'];

%% Create the cached file.
mkdir(dirName);
save(fileName,...
    'cbs_populationTotal'...
);

%% Remove clutter from the workspace
clear('dateName', 'fileName', 'dirName');