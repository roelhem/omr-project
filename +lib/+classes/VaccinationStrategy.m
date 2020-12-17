classdef VaccinationStrategy < lib.classes.AgeGroupPopulation
    %VACCINATIONSTRATEGY Holds the information of a vaccination stategy.
    %   Contains some methods to check if the VaccinationStrategy is
    %   allowed by the restrictions and methods to calculate the values
    %   that can be inserted into the SIRV_model.
    
    %% Initialisation
    methods
        function obj = VaccinationStrategy(scale, Rho, T)
            %VACCINATIONSTRATEGY Creates a new instance of this vaccination
            %   strategy.
            assert(isa(scale, 'lib.classes.AgeGroupPopulation'));
            obj = obj@lib.classes.AgeGroupPopulation(...
                scale.Boundaries, ...
                scale.Population ...
            );
            
            % Getting the default data.
            if nargin < 2
                Rho = zeros(obj.m, 0);
            end
            if nargin < 3
                T = 1;
            end
            
            if ~isduration(T)
                T = days(T);
            end
            
            % Initialize the timetable
            if isa(Rho, 'timetable')
                assert(width(Rho) == obj.m);
                obj.TTable = Rho;
            else
                assert(height(Rho) == obj.m);
                if width(Rho) == size(T)
                    obj.TTable = timetable( ...
                        'Size', [width(Rho) obj.m], ...
                        'VariableTypes', repmat("double", obj.m, 1), ...
                        'RowTimes', T ...
                    );
                elseif size(T) == 1
                    obj.TTable = timetable( ...
                        'Size', [width(Rho) obj.m], ...
                        'VariableTypes', repmat("double", obj.m, 1), ...
                        'TimeStep', T ...
                    );
                else
                    error('Invalid Time vector.');
                end
                
                obj.TTable{:,:} = Rho';
            end
            
            % Set the timetable variable names.
            obj.TTable.Properties.VariableNames = obj.GroupName;
            
            % Set the default InitialV.
            obj.InitialV = zeros(obj.m, 1);
        end
    end
    
    %% Stored properties
    properties
        TTable timetable % The underlying timetable holding the Vaccination
                         % Strategy data.
        InitialV double  % [m x 1] The initial amount of vaccinated people
                         %         per group j.
    end
    
    %% Dependent properties
    properties(Dependent)
        n      double    % [1 x 1] Number of timestamps for which a
                         %         value of the vaccination strategy is
                         %         known.
        Rho    double    % [m x n] Number of people that move from U to V
                         %         in group m at moment n due to this
                         %         vaccination strategy per day.
        RhoTot double    % [1 x n] The total of Rho over all groups.
        Nu     double    % [m x n] The nu that can be used in our model.
        P      double    % [m x n] The total amount of people that will
                         %         be vaccinated before the vaccination
                         %         moment.
        Ptot   double    % [m x 1] The total amount of people that will
                         %         be vaccinated at the and of this
                         %         vaccination strategy.
        T      duration  % [1 x n] The time after the initial state on
                         %         which the vaccination n will take
                         %         placed.
        DeltaT   double  % [1 x n] The amount of days from the previous
                         %         timestep, starting from 0.
        N        double  % [m x 1] Alias for the population of the age 
                         %         groups.
        InitialU double  % [m x 1] The initial amount of unvaccinated
                         %         people per group j.
        V        double  % [m x n] The amount of vaccinated people in
                         %         group j before datapoint n.
        U        double  % [m x n] The amount of unvaccinated people in
                         %         group j before datapoint n.
    end
    
    methods
        function out = get.n(obj)
            out = height(obj.TTable);
        end
        
        function out = get.Rho(obj)
            out = obj.TTable{:,:}';
        end
        
        function obj = set.Rho(obj, value)
            obj.TTable{:,:} = value';
        end
        
        function out = get.RhoTot(obj)
            out = sum(obj.Rho, 1);
        end
        
        function out = get.Nu(obj)
            out = obj.Rho ./ obj.U;
        end
        
        function out = get.P(obj)
            out = zeros(obj.m, obj.n);
            for i = 1:obj.n-1
                out(:,i+1) = out(:,i) + obj.Rho(:,i) .* obj.DeltaT(:,i);
            end
        end
        
        function obj = set.P(obj, value)
            assert(height(value) == obj.m);
            assert(width(value) == obj.n);
            for i = 1:obj.n-1
                obj.Rho(:,i) = (value(:,i+1) - value(:,i)) ./ obj.DeltaT(:,i);
            end
        end
        
        function out = get.Ptot(obj)
            out = sum(table2array(obj.DailyTable))';
        end
        
        function out = get.T(obj)
            out = obj.TTable.Properties.RowTimes';
        end
        
        function obj = set.T(obj, value)
            obj.TTable.Properties.RowTimes = value';
        end
        
        function out = get.DeltaT(obj)
            out = (obj.T(2:end) - obj.T(1:end-1)) / days(1);
        end
        
        function out = get.N(obj)
            out = obj.Population;
        end
        
        function obj = set.N(obj, value)
            assert(height(value) == obj.m);
            assert(width(value) == 1);
            obj.Population = value;
        end
        
        function out = get.InitialU(obj)
            out = obj.N - obj.InitialV;
        end
        
        function obj = set.InitialU(obj, value)
            assert(height(value) == obj.m);
            assert(width(value) == 1);
            obj.InitialV = obj.N - value;
        end
        
        function out = get.V(obj)
            out = repmat(obj.InitialV, 1, obj.n) + obj.P;
        end
        
        function obj = set.V(obj, value)
            obj.InitialV = value(:,1);
            obj.P = repmat(obj.InitialV, 1, obj.n) - value;
        end
        
        function out = get.U(obj)
            out = repmat(obj.N, 1, obj.n) - obj.V;
        end
        
        function obj = set.U(obj, value)
            obj.V = repmat(obj.N, 1, obj.n) - value;
        end
    end
    
    %% Reformatting the vaccination strategy.
    properties(Dependent)
        DailyTable duration % Timetable with a row for each day.
    end
    
    methods
        function out = get.DailyTable(obj)
            out = retime(obj.TTable, 'daily','mean');
            out = fillmissing(out, 'constant', 0);
        end
    end
    
    methods
        function result = V_at(obj, duration)
            % V_AT Returns the value of V at the provided duration.
            duration = days(duration);
            closestIndex = find(obj.T <= duration, 1, 'last');
            closestTime = obj.T(closestIndex);
            
            % Just give the InitialV if no times were available before
            % the provided time, just give the initial V.
            if(isempty(closestTime))
                result = obj.InitialV;
                return
            end
            
            % Compute the difference between the closest and the wanted
            % timestep.
            diffClosest = (duration) - closestTime;
            
            % If the difference is bigger than one day, it means that
            % we won't vaccinate on the days inbetween. Therefore, we can
            % just collapse the difference to one day.
            if(diffClosest >= days(1))
                diffClosest = days(1);
            end
            
            % Calculate the result.
            closestRho = obj.Rho(:,closestIndex);
            closestV   = obj.V(:,closestIndex);
            result = closestV + closestRho * diffClosest / days(1);
        end
        
        function result = retime(obj, NewT)
            % Get the new length of datapoints.
            if isduration(NewT)
                NewT = NewT/days(1);
            end
            New_n = length(NewT);
            
            % Compute the new Rho.
            NewRho = zeros(obj.m, New_n);
            for i = 1:New_n
                d = days(floor(NewT(i)));
                v = obj.DailyTable{d, :}';
                if ~isempty(v)
                    NewRho(:,i) = v;
                end
            end
            
            % Compute the new initial V.
            NewInitialV = obj.V_at(NewT(1));
            
            % Construct the new result.
            result = obj;
            result.InitialV = NewInitialV;
            result.TTable = array2timetable(NewRho', ...
                'RowTimes', days(NewT), ...
                'VariableNames', obj.GroupName ...
            );
        end
    end
    
end

