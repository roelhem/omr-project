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
    
    %% Checking Vaccination strategies.
    methods
        function out = checkPopLimit(obj, Strategy)
            out = all(Strategy.Ptot <= Strategy.N * obj.Efficacy);
        end
        
        function out = checkTotalMax(obj, Strategy)
            out = sum(Strategy.Ptot) <= obj.TotalMax * obj.Efficacy;
        end
        
        function out = checkMaxPerDay(obj, Strategy)
            out = all(Strategy.RhoTot <= obj.MaxPerDay * obj.Efficacy);
        end
        
        function out = check(obj, Strategy)
            out = obj.checkPopLimit(Strategy) && ...
                obj.checkTotalMax(Strategy) && ...
                obj.checkMaxPerDay(Strategy);
        end
    end
end

