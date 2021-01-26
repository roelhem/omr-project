function [InitialState] = SIRV_initial(varargin)
%SIRV_initial Creates a new initial state for the model.

% Getting the initialisation properties.
StartDate = datetime() - calmonths(1);
AgeGroups = lib.classes.ModelState.STD_Boundaries;
ContactMatrixFile = 'Contact_matrix.csv';
for ii = 1:2:length(varargin)
    switch string(varargin{ii})
        case "StartDate"
            StartDate = varargin{ii + 1};
        case "AgeGroups"
            AgeGroups = varargin{ii + 1};
        case "ContactMatrixFile"
            ContactMatrixFile = varargin{ii + 1};
    end
end

% Initialize the ModelState.
InitialState = lib.classes.ModelState.standardData( ...
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
        case "V"
            InitialState.V = varargin{ii + 1};
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


end

