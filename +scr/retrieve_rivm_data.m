
%% Loading the data
disp('Loading the data from the RIVM API...');

disp('Loading infectious people...');
rivm_infectiousPeople = lib.loaders.rivm_infectiousPeople;

disp('Loading reproduction values (R_0)...');
rivm_reproduction = lib.loaders.rivm_reproduction;

disp('Loading sewage study...');
rivm_sewage = lib.loaders.rivm_sewage;

disp('Loading totals per municipality study...');
rivm_totalPerMunicipality = lib.loaders.rivm_totalPerMunicipality;

disp('Loading covid cases...');
rivm_cases = lib.loaders.rivm_cases;

%% Determine the filename
dateName = datestr(datetime,"yyyy-mm-ddTHH-MM-ss");
dirName = 'cache/rivm/';
fileName = [dirName dateName '.mat'];

%% Create the cached file.
mkdir(dirName);
save(fileName,...
    'rivm_infectiousPeople',...
    'rivm_reproduction',...
    'rivm_sewage',...
    'rivm_totalPerMunicipality',...
    'rivm_cases'...
);