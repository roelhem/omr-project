classdef VaccinationRestriction
    %VACCINATIONRESTRICTION Contains the restrictions on a certain vaccin.
    
    %% Initialisation
    
    properties
        Efficacy  double % [1 x 1] The efficacy of the vaccine.
        MaxPerDay double % [1 x 1] The amount of vaccinations available
                         %         in one day.
        TotalMax  double % [1 x 1] The total amount of vaccines that can
                         %         be given.
    end
    
    methods
        function obj = VaccinationRestriction(Efficacy, MaxPerDay, TotalMax)
            %VACCINATIONRESTRICTION Construct an instance of this class
            
            % Get the default values.
            if nargin < 2
                Efficacy = 1;
            end
            if nargin < 1
                MaxPerDay = inf;
            end
            if nargin < 3
                TotalMax = inf;
            end
            
            % Set the default values.
            obj.Efficacy = Efficacy;
            obj.MaxPerDay = MaxPerDay;
            obj.TotalMax = TotalMax;
        end
    end
    
    %% Dependent properties.
    properties(Dependent)
        EffMaxPerDay double
        EffTotalMax  double
    end
    
    methods
        function out = get.EffMaxPerDay(obj)
            out = obj.MaxPerDay * obj.Efficacy;
        end
        
        function out = get.EffTotalMax(obj)
            out = obj.TotalMax * obj.Efficacy;
        end
    end
    
    %% Checking Vaccination strategies.
    methods
        function out = checkPopLimit(obj, Strategy)
            out = all(Strategy.Ptot <= Strategy.N * obj.Efficacy);
        end
        
        function out = checkTotalMax(obj, Strategy)
            out = sum(Strategy.Ptot) <= obj.EffTotalMax;
        end
        
        function out = checkMaxPerDay(obj, Strategy)
            out = all(Strategy.RhoTot <= obj.EffMaxPerDay);
        end
        
        function out = check(obj, Strategy)
            out = obj.checkPopLimit(Strategy) && ...
                obj.checkTotalMax(Strategy) && ...
                obj.checkMaxPerDay(Strategy);
        end
    end
    
    %% Strategy generator.
    methods
        function out = generateOrderStrategy(obj, Groups, order)
            global cbs_AgeGroupPopulation
            
            % Get the groups.
            if(~isa(Groups, 'lib.classes.AgeGroupPopulation'))
                Groups = cbs_AgeGroupPopulation.rescale(Groups);
            end
            
            if nargin < 3
                order = flip(1:Groups.m);
            end
            
            % Get some kind of upper bound for the init of out.
            maxVacc = min(Groups.PopulationTotal * obj.Efficacy, obj.EffTotalMax);
            
            % Alias the population, for better readability.
            N = Groups.Population;
            disp(maxVacc);
            disp(N);
            
            % Just vaccinate the maximum amount immediately when there is
            % no daily bound.
            if isinf(obj.MaxPerDay)
                out = zeros(Groups.m,1);
                bTot = 0;
                for i = order
                    disp(i);
                    a = min(maxVacc - bTot, N(i) * obj.Efficacy);
                    out(i) = a;
                    bTot = bTot + a;
                end
                return
            end
            
            % Get the maximum amount of days that the vaccination plan
            % might take, for initializing the result.
            maxDays = ceil(maxVacc / obj.EffMaxPerDay);
            
            % Initialize the results.
            out = zeros(Groups.m, maxDays);
            
            % Fill the values, vaccinating the people in the order first.
            bTot = 0;
            maxDay = obj.EffMaxPerDay;
            for j = 1:maxDays
                dTot = 0;
                for i = order
                    a = min([maxDay - dTot, maxVacc - bTot, N(i) * obj.Efficacy - sum(out(i,:))]);
                    out(i, j) = a;
                    dTot = dTot + a;
                    bTot = bTot + a;
                end
            end
            
        end
    end
end

