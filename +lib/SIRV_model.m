function M = SIRV_model(DeltaT, varargin)
%SIRV_MODEL Creates and runs a new model from standard data.
%   Detailed explanation goes here

% Set default arguments.
if nargin < 1
    DeltaT = 1;
end

% Get the optional arguments.
StartDate = datetime() - calmonths(1);
Method = "EulerForward";
n = 100;
AgeGroups = lib.classes.ModelState.STD_Boundaries;
ContactMatrixFile = 'Contact_matrix.csv';
for ii = 1:2:length(varargin)
    switch string(varargin{ii})
        case "StartDate"
            StartDate = varargin{ii + 1};
        case "Method"
            Method = varargin{ii + 1};
        case "Steps"
            n = varargin{ii + 1};
        case "n"
            n = varargin{ii + 1};
        case "AgeGroups"
            AgeGroups = varargin{ii + 1};
        case "ContactMatrixFile"
            ContactMatrixFile = varargin{ii + 1};
    end
end


% Get the initial state.
InitialState  = lib.classes.ModelState.standardData( ...
    StartDate, ...
    AgeGroups, ...
    ContactMatrixFile ...
);

% Modifying the initial state.
for ii = 1:2:length(varargin)
    switch string(varargin{ii})
        case "SympRatio"
            InitialState.SympRatio = varargin{ii + 1};
        case "AsympRatio"
            InitialState.AsympRatio = varargin{ii + 1};
        case "S"
            InitialState.S = varargin{ii + 1};
        case "I"
            InitialState.I = varargin{ii + 1};
        case "R"
            InitialState.R = varargin{ii + 1};
        case "Nu"
            InitialState.Nu = varargin{ii + 1};
        case "Alpha"
            InitialState.Alpha = varargin{ii + 1};
        case "Tau"
            InitialState.Tau = varargin{ii + 1};
        case "ReprNum"
            InitialState.ReprNum = varargin{ii + 1};
        case "ReprEff"
            InitialState.ReprEff = varargin{ii + 1};
    end
end

% Run the model to get an undated result.
UndatedResult = InitialState.run(DeltaT, n, Method);

% Get the result with dates.
M = UndatedResult.toDatedModelResult(StartDate);

end

