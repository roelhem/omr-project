classdef DatedModelResult < lib.classes.ModelResult
    %DATEDMODELRESULT A ModelResult that knows from which date the result
    %   was running.
    
    %% Initialisation
    properties
        StartDate datetime % [1 x 1] The date from which the model will 
                           %         start.
    end
    
    methods
        function obj = DatedModelResult(State, DeltaT, Method, StartDate)
            %DATEDMODELRESULT Creates a new Dated Model Result.
            obj = obj@lib.classes.ModelResult(State, DeltaT, Method);
            obj.StartDate = StartDate;
        end
    end
    
    %% Dependent properties
    properties(Dependent)
        Date datetime % [1 x n] The datetime that belongs to timestamp j.
    end
    
    methods
        function out = get.Date(obj)
            out = obj.StartDate + obj.T;
        end
    end
    
    %% Superclass overrides.
    methods
        function out = toTable(obj, cols, varargin)
            out = toTable@lib.classes.ModelResult(obj, cols,...
                'StartDate', obj.StartDate, ...
                varargin{:} ...
            );
        end
    end
end

