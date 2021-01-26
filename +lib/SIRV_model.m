function M = SIRV_model(DeltaT, varargin)
%SIRV_MODEL Creates and runs a new model from standard data.
%   Detailed explanation goes here

% Set default arguments.
if nargin < 1
    DeltaT = 1;
end

% Get the optional arguments.
Method = "EulerForward";
n = 100;
StartDate = datetime() - calmonths(1);
VaccinationStrategy = [];
VaccinationStrategyTimeStep = 1;
InitialState = [];
for ii = 1:2:length(varargin)
    switch string(varargin{ii})
        case "Method"
            Method = varargin{ii + 1};
        case "Steps"
            n = varargin{ii + 1};
        case "StartDate"
            StartDate = varargin{ii + 1};
        case "n"
            n = varargin{ii + 1};
        case "VaccinationStrategy"
            VaccinationStrategy = varargin{ii + 1};
        case "VaccinationStrategyTimeStep"
            VaccinationStrategyTimeStep = varargin{ii + 1};
        case "InitialState"
            InitialState = varargin{ii + 1};
    end
end


% Get the default initial state
if isempty(InitialState)
    InitialState = lib.SIRV_initial(varargin{:});
end

% Get the vaccination strategy.]
if isempty(VaccinationStrategy)
    VaccinationStrategy = InitialState.toVaccinationStrategy();
elseif isa(VaccinationStrategy, 'lib.classes.VaccinationStrategy')
    VaccinationStrategy.InitialV = InitialState.V;
elseif isnumeric(VaccinationStrategy)
    VaccinationStrategy = lib.classes.VaccinationStrategy( ...
        InitialState, ...
        VaccinationStrategy, ...
        VaccinationStrategyTimeStep ...
    );
end

% Run the model to get an undated result.
UndatedResult = InitialState.run(DeltaT, n, ...
    Method, ...
    VaccinationStrategy ...
);

% Get the result with dates.
M = UndatedResult.toDatedModelResult(StartDate);

end

