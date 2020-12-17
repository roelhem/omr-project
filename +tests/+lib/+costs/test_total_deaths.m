% Load necessary data into workspace
scr.init
%% Get cost parameters

% Input lethality vector
l = [0.00022274,0.00022274,0.00022274,0.00022274,0.00022274,0.00221398,0.01351613,0.134212];

% Run model
% Call the lib.SIRV_model.
M = lib.SIRV_model();

%% Calculate total deaths
lib.costs.total_deaths(M,l)