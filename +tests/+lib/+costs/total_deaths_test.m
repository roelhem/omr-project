% Load necessary data into workspace
scr.init
%% Get cost parameters

% Input lethality vector
l = [0.00022274,0.00022274,0.00022274,0.00022274,0.00022274,0.00221398,0.01351613,0.134212];

% Run model
% Call the lib.SIRV_model.
M = lib.SIRV_model(1, ... The timestap DeltaT
    'StartDate', '2020-10-09', ... The date for which the data should be loaded.
    'Steps', 160 , ... The amount of steps to calculate beforehand.
    'Method', "EulerForward", ... The numeric method to solve the ODE of the model.
    'AgeGroups', 0:10:70, ... The way the age groups should be devided.
    'ContactMatrixFile', "Contact_matrix.csv", ... The file that contains the contact matrix.
    'Tau', 5, ... The average time of an infection.
    'ReprNum', 5.0 ... Set the reproduction number to a fixed value.
);

%% Calculate total deaths
lib.costs.total_deaths(M,l)