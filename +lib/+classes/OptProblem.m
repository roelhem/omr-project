classdef OptProblem
    %OPTPROBLEM Holds all the parameters for an optimisation problem.
    
    
    %% Initialisation
    properties
        InitialState           lib.classes.ModelState
        VaccinationRestriction lib.classes.VaccinationRestriction
        CostFunc
        DeltaT                 double
        Method                 string
        Steps                  double
    end
    
    methods
        function obj = OptProblem(InitialState, VaccinationRestriction, ...
               CostFunc, DeltaT, Steps, Method)
            %OPTPROBLEM Construct an instance of this class
            
            if nargin < 3
                CostFunc = @lib.costs.total_deaths;
            end
            if nargin < 4
                DeltaT = 1;
            end
            if nargin < 5
                Steps = ceil(VaccinationRestriction.expectedDays(InitialState) / DeltaT);
            end
            if nargin < 6
                Method = 'EulerForward';
            end
            
            obj.VaccinationRestriction = VaccinationRestriction;
            obj.InitialState = InitialState;
            obj.CostFunc = CostFunc;
            obj.DeltaT = DeltaT;
            obj.Steps = Steps;
            obj.Method = Method;
        end
    end
    
    %% Dependend properties
    properties
        m double
        days double
        l double
    end
    
    methods
        function out = get.m(obj)
            out = obj.InitialState.m;
        end
        
        function out = get.days(obj)
            out = floor(obj.Steps * obj.DeltaT);
        end
        
        function out = get.l(obj)
            out = obj.m * obj.days;
        end
    end
    
    %% Restrictions
    methods
        function [A,b] = solutionSpace(obj)
            % Population restriction.
            A_pop = zeros(obj.m, obj.l);
            b_pop = zeros(obj.m, 1);
            for i = 1:obj.m
                A_pop(i, i:obj.m:obj.l) = 1;
                b_pop(i) = obj.InitialState.N(i) * obj.VaccinationRestriction.Efficacy;
            end
            
            % Start building the result.
            A = A_pop;
            b = b_pop;
            
            % Vaccinations per day.
            if ~isinf(obj.VaccinationRestriction.EffMaxPerDay)
                A_daily = zeros(obj.days, obj.l);
                b_daily = ones(obj.days,1) * obj.VaccinationRestriction.EffMaxPerDay;
                for i = 1:obj.days
                    startIndex = obj.m * (i - 1) + 1;
                    endIndex = obj.m * i;
                    A_daily(i, startIndex:endIndex) = 1;
                end
                A = [A;A_daily];
                b = [b;b_daily];
            end
            
            % Restriction on total
            if ~isinf(obj.VaccinationRestriction.EffTotalMax)
                A_total = ones(1, obj.l);
                b_total = obj.VaccinationRestriction.EffTotalMax;
                A = [A;A_total];
                b = [b;b_total];
            end
            
            % Positive values restriction.
            A_pos = -eye(obj.l, obj.l);
            b_pos = zeros(obj.l,1);
            A = [A;A_pos];
            b = [b;b_pos];
            
        end
    end
    
    %% Helper functions
    methods
        function M = input2stratmat(obj, X)
            assert(length(X) == obj.l);
            M = reshape(X,obj.m,obj.days);
        end
        
        function X = stratmat2input(obj, M)
            X = reshape(M,obj.l,1);
        end
    end
    
    %% Optimisation function
    methods
        function out = getOptimisationHandle(obj, l)
            
            out = @(x) obj.CostFunc(obj.InitialState.run( ...
                obj.DeltaT, ...
                obj.Steps, ...
                obj.Method, ...
                lib.classes.VaccinationStrategy( ...
                    obj.InitialState, ...
                    obj.input2stratmat(x) ...
                ) ...
            ), l);
            
        end
    end
end

