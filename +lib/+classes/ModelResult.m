classdef ModelResult < lib.classes.AgeGroupPopulation
    %MODELRESULT The results got form a model that run successfully.
    
    %% Initialisation
    methods
        function obj = ModelResult(State, DeltaT, Method)
            %MODELRESULT Construct an instance of this class.
            
            % Check if there were states given.
            if isempty(State)
                error("Empty result given");
            end
            
            % Initialize the object.
            obj = obj@lib.classes.AgeGroupPopulation(...
                State(1).Boundaries, ...
                State(1).Population ...
            );
            
            % Initialize the properties.
            obj.State = State;
            obj.DeltaT = DeltaT;
            obj.Method = categorical(Method);
        end
    end
    
    methods(Static)
        function obj = fromInitialState(InitialState, DeltaT, n, Method)
            % FROMINITIAL Creates a new ModelResult by running the model
            %             starting from an initial state.
            
            % Checking/parsing the input.
            assert(isa(InitialState, 'lib.classes.ModelState'), 'InitialState have to be a lib.classes.ModelState');
            obj = InitialState.run(DeltaT, n, Method);
        end
    end
    
    %% Stored properties.
    properties
        State    lib.classes.ModelState % [1 x n] Array of model states.
        DeltaT   double                 % [1 x 1] Size of each timestep.
        Method   categorical            % Method that was used to compute
                                        % the result of the model.
    end
    
    %% Dependent properties.
    properties(Dependent)
        n double   % [1 x 1] The amount of steps in this model.
    end
    
    methods
        function out = get.n(obj)
            out = width(obj.State);
        end
        
        function obj = set.n(obj, value)
            obj = obj.resize(value);
        end
    end
    
    %% Conversion to other data types.
    methods
        function out = toTable(obj, cols, varargin)
            columns = struct();
            
            % Getting the optional parameters
            GroupSubTables = false;
            TimeTable = true;
            StartDate = [];
            for ii = 1:2:length(varargin)
                switch string(varargin{ii})
                    case "GroupSubTables"
                        GroupSubTables = varargin{ii + 1};
                    case "TimeTable"
                        TimeTable = varargin{ii + 1};
                    case "StartDate"
                        StartDate = varargin{ii + 1};
                end
            end
            
            % Contstructing the table.
            for i = 1:length(cols)
                colname = cols{i};
                switch string(colname)
                    case "T"
                        columns.T = obj.T';
                    otherwise
                        value = obj.getValue(colname)';
                        if GroupSubTables && width(value) == obj.m
                            subTable = array2table(value);
                            subTable.Properties.VariableNames = obj.GroupName;
                            columns.(colname) = subTable;
                        else
                            columns.(colname) = value;
                        end
                end
            end
            
            % Get the table.
            out = struct2table(columns);
            if TimeTable
                if ~isempty(StartDate)
                    out = table2timetable(out, 'RowTimes', datetime(StartDate) + obj.T');
                else
                    out = table2timetable(out, 'RowTimes', obj.T');
                end
            end
        end
        
        function result = toDatedModelResult(obj, StartDate)
            result = lib.classes.DatedModelResult( ...
                obj.State, ...
                obj.DeltaT, ...
                obj.Method, ...
                StartDate ...
            );
        end
    end
    
    %% Derrived properties from the ModelState.
    properties(Dependent)
        T duration % [1 x n] The duration of the timestep from he initail
                   %         state.
        N          double % [m x n]
        S          double % [m x n]
        I          double % [m x n]
        R          double % [m x n]
        U          double % [m x n]
        V          double % [m x n]
        Alpha      double % [m x n]
        Tau        double % [m x n]
        Beta       double % [m x m x n]
        Nu         double % [m x n]
        SusVect    double % [m x n]
        C          double % [m x m x n]
        P_CtoI     double % [m x m x n]
        ReprNum    double % [1 x n]
        ReprEff    double % [1 x n]
        SympRatio  double % [m x n]
        AsympRatio double % [m x n]
        I_Symp     double % [m x n]
        I_Asymp    double % [m x n]
        StoI       double % [m x n]
        ItoR       double % [m x n]
        StoV       double % [m x n]
        ItoV       double % [m x n]
        RtoV       double % [m x n]
        UtoV       double % [m x n]
        dS         double % [m x n]
        dI         double % [m x n]
        dR         double % [m x n]
        dU         double % [m x n]
        dV         double % [m x n]
    end
    
    methods
        function out = get.T(obj)
            out = obj.getValue('T');
        end
        
        function out = get.N(obj)
            out = obj.getValue('N');
        end
        
        function out = get.S(obj)
            out = obj.getValue('S');
        end
        
        function out = get.I(obj)
            out = obj.getValue('I');
        end
        
        function out = get.R(obj)
            out = obj.getValue('R');
        end
        
        function out = get.U(obj)
            out = obj.getValue('U');
        end
        
        function out = get.V(obj)
            out = obj.getValue('V');
        end
        
        function out = get.Alpha(obj)
            out = obj.getValue('Alpha');
        end
        
        function out = get.Tau(obj)
            out = obj.getValue('Alpha');
        end
        
        function out = get.Beta(obj)
            out = obj.getValue('Beta');
        end
        
        function out = get.Nu(obj)
            out = obj.getValue('Nu');
        end
        
        function out = get.SusVect(obj)
            out = obj.getValue('SusVect');
        end
        
        function out = get.C(obj)
            out = obj.getValue('C');
        end
        
        function out = get.P_CtoI(obj)
            out = obj.getValue('P_CtoI');
        end
        
        function out = get.ReprNum(obj)
            out = obj.getValue('ReprNum');
        end
        
        function out = get.ReprEff(obj)
            out = obj.getValue('ReprEff');
        end
        
        function out = get.SympRatio(obj)
            out = obj.getValue('SympRatio');
        end
        
        function out = get.AsympRatio(obj)
            out = obj.getValue('AsympRatio');
        end
        
        function out = get.I_Symp(obj)
            out = obj.getValue('I_Symp');
        end
        
        function out = get.I_Asymp(obj)
            out = obj.getValue('I_Asymp');
        end
        
        function out = get.StoI(obj)
            out = obj.getValue('StoI');
        end
        
        function out = get.ItoR(obj)
            out = obj.getValue('ItoR');
        end
        
        function out = get.StoV(obj)
            out = obj.getValue('StoV');
        end
        
        function out = get.ItoV(obj)
            out = obj.getValue('StoV');
        end
        
        function out = get.RtoV(obj)
            out = obj.getValue('RtoV');
        end
        
        function out = get.UtoV(obj)
            out = obj.getValue('UtoV');
        end
    end
    
    %% Helper methods.
    methods
        function out = getValue(obj, fieldName)
            c = arrayfun(@(x)x.(fieldName), obj.State, 'UniformOutput', false);
            if(width(c{1}) == 1)
                % Init the result.
                out = zeros(height(c{1}), obj.n);
                % Check if the elements are of type duration.
                if isduration(c{1})
                    out = days(out);
                end
                for i = 1:obj.n
                    out(:,i) = c{i};
                end
            else
                out = zeros(width(c{1}), height(c{1}), obj.n);
                for i = 1:obj.n
                    out(:,:,i) = c{i};
                end
            end
        end
    end
    
    %% Extending the result.
    methods
        function result = extend(obj, n)
            newResult = obj.State(end).run(obj.DeltaT, n + 1, obj.Method);
            result = obj;
            result.State(end+1:end+n) = newResult.State(2:end);
        end
        
        function result = shrink(obj, n)
            result = obj;
            result.State = result.State(1:end-n);
        end
        
        function result = resize(obj, n)
            dn = n - obj.n;
            if dn > 0
                result = obj.extend(dn);
            elseif dn < 0
                result = obj.shrink(-dn);
            else
                result = obj;
            end
        end
    end
end

