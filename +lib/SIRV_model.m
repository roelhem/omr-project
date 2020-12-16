function M = SIRV_model(DeltaT, n, varargin)
%SIRV_MODEL Creates and runs a new model from standard data.
%   Detailed explanation goes here

% Set default arguments.
if nargin < 1
    DeltaT = 1;
end
if nargin < 2
    n = 100;
end

StartDate = datetime() - calmonths(1);
Method = "EulerForward";
for ii = 1:2:length(varargin)
    switch string(varargin{ii})
        case "StartDate"
            StartDate = varargin{ii + 1};
        case "Method"
            Method = varargin{ii + 1};
    end
end



InitialState  = lib.classes.ModelState.standardData(StartDate);
UndatedResult = InitialState.run(DeltaT, n, Method);
M = UndatedResult.toDatedModelResult(StartDate);

end

