classdef VaccinationRestriction < lib.classes.Serializable
    %VACCINATIONRESTRICTION Contains the restrictions on a certain vaccin.
    
    %% Serialisation
    methods
        function out = toStruct(obj)
            out = struct(...
                'VaccineEfficacy', obj.Efficacy, ...
                'MaxVaccinesPerDay', obj.MaxPerDay, ...
                'MaxVaccinesTotal', obj.TotalMax, ...
                'ConstraintTolerance', obj.ConstraintTolerance ...
            );
        end
    end
    
    methods(Static)
        function out = fromStruct(s)
            out = lib.classes.VaccinationRestriction(...
                'VaccineEfficacy', s.VaccineEfficacy, ...
                'MaxVaccinesPerDay', s.MaxVaccinesPerDay, ...
                'MaxVaccinesTotal', s.MaxVaccinesTotal, ...
                'Tolerance', s.ConstraintTolerance ...
            );
        end
    end
    
    %% Initialisation
    
    properties
        Efficacy  double % [1 x 1] The efficacy of the vaccine.
        MaxPerDay double % [1 x 1] The amount of vaccinations available
                         %         in one day.
        TotalMax  double % [1 x 1] The total amount of vaccines that can
                         %         be given.
        ConstraintTolerance double % the tolerance of the boundaries.
    end
    
    methods
        function obj = VaccinationRestriction(varargin)
            %VACCINATIONRESTRICTION Construct an instance of this class
            
            % Get the default values.
            obj.Efficacy = 1;
            obj.MaxPerDay = inf;
            obj.TotalMax = inf;
            obj.ConstraintTolerance = 1e-3;
            
            % Set the values from the arguments.
            for ii = 1:2:length(varargin)
                switch string(varargin{ii})
                    case "VaccineEfficacy"
                        obj.Efficacy = varargin{ii + 1};
                    case "MaxVaccinesPerDay"
                        obj.MaxPerDay = varargin{ii + 1};
                    case "MaxVaccinesTotal"
                        obj.TotalMax = varargin{ii + 1};
                    case "Tolerance"
                        obj.ConstraintTolerance = varargin{ii + 1};
                end
            end

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
            out = all(Strategy.Ptot <= Strategy.N * obj.Efficacy + obj.ConstraintTolerance);
        end
        
        function out = checkTotalMax(obj, Strategy)
            out = sum(Strategy.Ptot) <= obj.EffTotalMax + obj.ConstraintTolerance;
        end
        
        function out = checkMaxPerDay(obj, Strategy)
            out = all(Strategy.RhoTot <= obj.EffMaxPerDay + obj.ConstraintTolerance);
        end
        
        function out = check(obj, Strategy)
            out = obj.checkPopLimit(Strategy) && ...
                obj.checkTotalMax(Strategy) && ...
                obj.checkMaxPerDay(Strategy);
        end
    end
    
    %% Helper methods
    methods
        function out = vaccUpperBound(obj, Groups)
            % Get some kind of upper bound for the init of out.
            out = min(Groups.PopulationTotal * obj.Efficacy, obj.EffTotalMax);
        end
        function out = expectedDays(obj, Groups)
            % Get some kind of upper bound for the init of out.
            out = ceil(obj.vaccUpperBound(Groups) / obj.EffMaxPerDay);
        end
    end
    
    methods
        
        function [vect, last] = getFillVect(obj, nvacc, initLim)
            maxDay = obj.EffMaxPerDay;
            first_val = min(nvacc, maxDay - initLim);
            next_nvacc = nvacc - first_val;
            full_days = floor(next_nvacc / maxDay);
            last = next_nvacc - full_days * maxDay;
            vect = [first_val, repmat(maxDay, 1, full_days), last];
        end
    end
    
    %% Strategy generator.
    methods
        
        function out = generateFullStrategyMat(obj, Groups, n)
            maxPerDay = obj.EffMaxPerDay;
            PE = Groups.Population * obj.Efficacy;
            
            % Initialise output.
            out = zeros(Groups.m, n);
            
            % Fill values.
            for i = 1:n
                % Amount of people to (effectively) vaccinate.
                UE = PE - sum(out, 2);
                if all(UE == 0)
                    break
                end
                
                % First distribution
                while sum(out(:,i)) < maxPerDay && any(UE > out(:, i))
                    % Compute the amount of vaccines that need to be
                    % distributed.
                    underflow = maxPerDay - sum(out(:,i));
                    
                    % Determine how to distribute them.
                    dist = rand(Groups.m, 1) .* max(0, UE - out(:,i));
                    dist = dist ./ sum(dist);
                    
                    % Compute optimistic distribution.
                    delta = round(dist * underflow);
                    err_size = underflow - sum(delta);
                    err_index = randi([1 Groups.m]);
                    delta(err_index) = delta(err_index) + err_size;
                    
                    % add to distribution.
                    out(:,i) = out(:,i) + min(UE - out(:,i), delta);
                end
            end
        end
        
        function out = generateFullStrategy(obj, Groups, n)
            global cbs_AgeGroupPopulation
            
            % Get the groups.
            if(~isa(Groups, 'lib.classes.AgeGroupPopulation'))
                Groups = cbs_AgeGroupPopulation.rescale(Groups);
            end
            
            if nargin < 3
                n = obj.expectedDays(Groups);
            end
            
            out = obj.generateFullStrategyMat(Groups, n);
            
            % Create the vaccination strategy as the result.
            out = lib.classes.VaccinationStrategy(Groups, out);
        end
        
        function out = generateUniformStrategyMat(obj, Groups, n)
            N = Groups.Population;
            N_eff = N * obj.Efficacy;
            
            % Init the result.
            out = zeros(Groups.m,n);
            
            % Just vaccinate the maximum amount immediately when there is
            % no daily bound.
            if isinf(obj.MaxPerDay)
                out(:,1) = N_eff;
                return
            end
            
            eta = obj.EffMaxPerDay / sum(N_eff);
            fullDays = min(n,floor(1 / eta));
            out(:,1:fullDays) = repmat(N_eff * eta, 1, fullDays);
            if fullDays < n
                out(:,fullDays + 1) = N_eff - sum(out, 2);
            end
        end
        
        function out = generateUniformStrategy(obj, Groups, n)
            global cbs_AgeGroupPopulation
            
            % Get the groups.
            if(~isa(Groups, 'lib.classes.AgeGroupPopulation'))
                Groups = cbs_AgeGroupPopulation.rescale(Groups);
            end
            
            out = obj.generateUniformStrategyMat(Groups, n);
            
            % Create the VaccinationStrategy object and return the result.
            out = lib.classes.VaccinationStrategy(Groups, out);
            
        end
        
        function out = generateOrderStrategyMat(obj, Groups, order, n)
            % Get some kind of upper bound for the init of out.
            maxVacc = obj.vaccUpperBound(Groups);
            
            % Alias the population, for better readability.
            N = Groups.Population;
            
            % Just vaccinate the maximum amount immediately when there is
            % no daily bound.
            if isinf(obj.MaxPerDay)
                out = zeros(Groups.m,n);
                bTot = 0;
                for i = order
                    a = min(maxVacc - bTot, N(i) * obj.Efficacy);
                    out(i) = a;
                    bTot = bTot + a;
                end
                return
            end
            
            % Get the maximum amount of days that the vaccination plan
            % might take, for initializing the result.
            if nargin < 4
                n = obj.expectedDays(Groups);
            end
            % Initialize the results.
            out = zeros(Groups.m, n);
            
            % Fill the values, vaccinating the people in the order first.
            bTot = 0;
            initLim = 0;
            i = 1;
            for j = order
                nvacc = min(maxVacc - bTot, N(j) * obj.Efficacy);
                if nvacc <= 0
                    break
                end
                [vect, last] = obj.getFillVect(nvacc, initLim);
                l = min(length(vect), n-i+1);
                out(j, i:i+l-1) = vect(1:l);
                i = i+l-1;
                bTot = bTot + nvacc;
                initLim = last;
            end
        end
        
        function out = generateOrderStrategy(obj, Groups, order, n)
            global cbs_AgeGroupPopulation
            
            % Get the groups.
            if(~isa(Groups, 'lib.classes.AgeGroupPopulation'))
                Groups = cbs_AgeGroupPopulation.rescale(Groups);
            end
            
            if nargin < 3
                order = flip(1:Groups.m);
            end
            
            out = obj.generateOrderStrategyMat(Groups, order, n);
            
            % Create the VaccinationStrategy object and return the result.
            out = lib.classes.VaccinationStrategy(Groups, out);
            
        end
       
    end
end

