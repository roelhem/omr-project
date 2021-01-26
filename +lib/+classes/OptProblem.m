classdef OptProblem < lib.classes.Serializable
    %OPTPROBLEM Holds all the parameters for an optimisation problem.
    
    %% Serialisation
    methods
        function out = toStruct(obj)
            out = struct(...
                'InitialState', obj.InitialState.toStruct(), ...
                'VaccinationRestriction', obj.VaccinationRestriction.toStruct(), ...
                'CostFunc', obj.CostFunc, ...
                'CostFuncParams', {obj.CostFuncParams}, ...
                'DeltaT', obj.DeltaT, ...
                'Method', obj.Method, ...
                'Steps', obj.Steps, ...
                'StartDate', obj.StartDate ...
            );
        end
    end
    
    methods(Static)
        function out = fromStruct(s)
            out = lib.classes.OptProblem(...
                'InitialState', lib.classes.ModelState.fromStruct(s.InitialState), ...
                'VaccinationRestriction', lib.classes.VaccinationRestriction.fromStruct(s.VaccinationRestriction), ...
                'CostFunc', s.CostFunc, ...
                'CostFuncParams', s.CostFuncParams, ...
                'DeltaT', s.DeltaT, ...
                'Method', s.Method, ...
                'Steps', s.Steps, ...
                'StartDate', s.StartDate ...
            );
        end
    end
    
    %% Initialisation
    properties
        InitialState           lib.classes.ModelState
        VaccinationRestriction lib.classes.VaccinationRestriction
        CostFunc
        CostFuncParams         cell
        DeltaT                 double
        Method                 string
        Steps                  double
        StartDate              datetime
    end
    
    methods
        function obj = OptProblem(varargin)
            %OPTPROBLEM Construct an instance of this class
            obj.VaccinationRestriction = lib.classes.VaccinationRestriction(varargin{:});
            obj.InitialState = lib.SIRV_initial(varargin{:});
            obj.CostFunc = 'total_deaths';
            obj.CostFuncParams = {};
            obj.DeltaT = 1;
            obj.Steps = [];
            obj.Method = 'EulerForward';
            obj.StartDate = datetime('2021-01-01');
            
            % Setting variables from varargin.
            for ii = 1:2:length(varargin)
                switch string(varargin{ii})
                    case "Method"
                        obj.Method = varargin{ii + 1};
                    case "Steps"
                        obj.Steps = varargin{ii + 1};
                    case "VaccinationRestriction"
                        obj.VaccinationRestriction = varargin{ii + 1};
                    case "CostFunc"
                        obj.CostFunc = varargin{ii + 1};
                    case "CostFuncParams"
                        if iscell(varargin{ii+1})
                            obj.CostFuncParams = varargin{ii+1};
                        else
                            obj.CostFuncParams = varargin(ii+1);
                        end
                    case "StartDate"
                        obj.StartDate = varargin{ii + 1};
                    case "InitialState"
                        obj.InitialState = varargin{ii + 1};
                end
            end
            
            
            % Default value for steps.
            if isempty(obj.Steps)
                obj.Steps = floor(obj.VaccinationRestriction.expectedDays(obj.InitialState) / obj.DeltaT) -1;
            end
        end
    end
    
    %% Dependend properties
    properties(Dependent)
        N double
        m double
        days double
        nvars double % Number of variables.
        Aineq_population double % Population restriction matrix.
        Bineq_population double % Population restriction boundary.
        Aineq_per_day double % Vaccinations per day restriction matrix.
        Bineq_per_day double % Vaccinations per day restriction boundary.
        Aineq_total double % Total vaccinations restriction matrix.
        Bineq_total double % Total vaccinations restriction boundary.
        Aineq_positive double % Positive restriction matrix.
        Bineq_positive double % Positive restriction boundary.
        Aineq double % Linear inequality restrictions matrix.
        Bineq double % Linear inequality restrictions boundary.
        Aeq double % Linear equality restrictions matrix.
        Beq double % Linear equality restrictions boundary.
        nonlcon % Non-linear constraints.
        lb double % Lower bound.
        ub double % Upper bound.
        fitnessfcn function_handle % The fitness function.
        options % The optimisation options.
        CreationFcn function_handle % Function that generates intial populations.
        MutationFcn function_handle % Function that add mutations to the population.
    end
    
    methods
        function out = get.N(obj)
            out = obj.InitialState.N;
        end
        
        function out = get.m(obj)
            out = obj.InitialState.m;
        end
        
        function out = get.days(obj)
            out = floor(obj.Steps * obj.DeltaT);
        end
        
        function out = get.nvars(obj)
            out = obj.m * obj.days;
        end
        
        function out = get.Aineq_population(obj)
            out = zeros(obj.m, obj.nvars);
            for i = 1:obj.m
                out(i, i:obj.m:obj.nvars) = 1;
            end
        end
        
        function out = get.Bineq_population(obj)
            out = zeros(obj.m, 1);
            for i = 1:obj.m
                out(i) = obj.InitialState.N(i) * obj.VaccinationRestriction.Efficacy;
            end
        end
        
        function out = get.Aineq_per_day(obj)
            if isinf(obj.VaccinationRestriction.EffMaxPerDay)
                % No per day restriction when EffMaxPerDay is infinite.
                out = zeros(0, obj.nvars);
            else
                out = zeros(obj.days, obj.nvars);
                for i = 1:obj.days
                    startIndex = obj.m * (i - 1) + 1;
                    endIndex = obj.m * i;
                    out(i, startIndex:endIndex) = 1;
                end
            end
        end
        
        function out = get.Bineq_per_day(obj)
            % No per day restriction when EffMaxPerDay is infinite.
            if isinf(obj.VaccinationRestriction.EffMaxPerDay)
                out = zeros(0, 1);
            else
                out = ones(obj.days,1) * obj.VaccinationRestriction.EffMaxPerDay;
            end
        end
        
        function out = get.Aineq_total(obj)
            % No per day restrictions when Total is infinite.
            if isinf(obj.VaccinationRestriction.EffTotalMax)
                out = zeros(0, obj.nvars);
            else
                out = ones(1, obj.nvars);
            end
        end
        
        function out = get.Bineq_total(obj)
            % No per day restrictions when Total is infinite.
            if isinf(obj.VaccinationRestriction.EffTotalMax)
                out = zeros(0, 1);
            else
                out = obj.VaccinationRestriction.EffTotalMax;
            end
        end
        
        function out = get.Aineq_positive(obj)
            out = -eye(obj.nvars, obj.nvars);
        end
        
        function out = get.Bineq_positive(obj)
            out = zeros(obj.nvars, 1);
        end
        
        function out = get.Aineq(obj)
            out = [obj.Aineq_population; obj.Aineq_per_day; obj.Aineq_total];
        end
        
        function out = get.Bineq(obj)
            out = [obj.Bineq_population; obj.Bineq_per_day; obj.Bineq_total];
        end
        
        function out = get.Aeq(obj)
            out = zeros(0, obj.nvars);
        end
        
        function out = get.Beq(obj)
            out = zeros(0, obj.nvars);
        end
        
        function out = get.nonlcon(obj)
            function [c, ceq] = f(x)
                c = [];
                ceq = [];
            end
            out = @f;
        end
        
        function out = get.lb(obj)
            out = zeros(1, obj.nvars);
        end
        
        function out = get.ub(obj)
            out = ones(obj.nvars, 1) * obj.VaccinationRestriction.EffMaxPerDay;
        end
        
        function out = get.fitnessfcn(obj)
            if string(obj.CostFunc) == "total_deaths"
                out = lib.costs.get_totalDeaths_vectorized(obj.InitialState, obj.CostFuncParams{:});
                return
            elseif string(obj.CostFunc) == "max_I"
                out = lib.costs.get_maxI_vectorized(obj.InitialState, obj.CostFuncParams{:});
                return
            elseif string(obj.CostFunc) == "Reff"
                out = lib.costs.get_Reff_vectorized(obj.InitialState, obj.CostFuncParams{:});
                return
            end
            
            vacc_fcn = obj.InitialState.vacc_vectorized_run;
            out = @(x)obj.CostFunc(vacc_fcn(x), obj.CostFuncParams{:});
        end
        
        function out = get.MutationFcn(obj)
            
            function MutatedChildren = mut_mean(parents, options, nvars, FitnessFcn, state, thisScore, thisPopulation)
                MutatedChildren = thisPopulation(parents,:);
                
                max_size = 20;
                
                for i = 1:length(parents)
                    child = obj.input2stratmat(thisPopulation(parents(i), :));
                    size = randi(max_size);
                    start = randi([1, obj.days - size]);
                    child(:, start:start+size-1) = repmat(mean(child(:, start:start+size-1), 2), 1, size);
                    MutatedChildren(i, :) = obj.stratmat2input(child);
                end
                
            end
            
            function MutatedChildren = mut_swap(parents, options, nvars, FitnessFcn, state, thisScore, thisPopulation)
                MutatedChildren = thisPopulation(parents,:);
        
                ps = height(parents);
                for i = 1:ps
                    parent = obj.input2stratmat(thisPopulation(parents(i), :));
                    available_days = sum(sum(parent, 1) > 1e-5);
                    
                    swap_size = randi([1 randi(6)]);
                    a = randi([swap_size+1, available_days - swap_size + 1]);
                    b = randi([1, a - swap_size]);
                    a = a:a+swap_size-1;
                    b = b:b+swap_size-1;
                    temp = parent(:, a);
                    parent(:, a) = parent(:, b);
                    parent(:, b) = temp;
                    MutatedChildren(i,:) = obj.stratmat2input(parent);
                    
%                     I_a = (1:obj.m) + (swap_vals(1) - 1)*obj.m;
%                     I_b = (1:obj.m) + (swap_vals(2) - 1)*obj.m;
%                     temp = MutatedChildren(psi, I_a);
%                     MutatedChildren(psi, I_a) = MutatedChildren(psi, I_b);
%                     MutatedChildren(psi, I_b) = temp;
                end
            end
            
            function MutatedChildren = mut_dist(parents, options, nvars, FitnessFcn, state, thisScore, thisPopulation)
                MutatedChildren = thisPopulation(parents,:);
                
                max_change = obj.VaccinationRestriction.EffMaxPerDay;
                
                for i = 1:length(parents)
                    parent = obj.input2stratmat(thisPopulation(parents(i), :));
                    
                    % Group to change
                    grp = randi(obj.m);
                    
                    % maximum change.
                    delta = unifrnd(0, max_change);
                    exp_days = obj.VaccinationRestriction.expectedDays(obj.InitialState);
                    
                    % Add distribution.
                    add_mu = unifrnd(1, exp_days);
                    add_sigma = unifrnd(0.2, 5);
                    add_dist = normpdf(1:obj.days, add_mu, add_sigma);
                    
                    % Compensate distribution.
                    cmp_mu = unifrnd(1, exp_days);
                    cmp_sigma = unifrnd(5, 10);
                    cmp_dist = normpdf(1:obj.days, cmp_mu, cmp_sigma);
                    
                    % Create mutation.
                    child = lib.mutations.mut_dist_base(parent, grp, delta, add_dist, cmp_dist);
                    MutatedChildren(i, :) = obj.stratmat2input(child);
                end
            end
            
            function MutatedChildren = f(parents, options, nvars, FitnessFcn, state, thisScore, thisPopulation)
                %randomStrategies = obj.getRandomFullStrategies(p);
                %gamma = diag(unifrnd(0, 0.5, p, 1));
                MutatedChildren = thisPopulation(parents,:);
                
                for i = 1:length(parents)
                    parent = thisPopulation(parents(i),:);
                    
                    % Get the value that will be mutated.
                    day_a = randi(obj.days);
                    day_b = randi(obj.days);
                    grp_a = randi(obj.m);
                    grp_b = randi(obj.m);
                    
                    aa_index = (day_a - 1)*obj.m + grp_a;
                    ab_index = (day_a - 1)*obj.m + grp_b;
                    ba_index = (day_b - 1)*obj.m + grp_a;
                    bb_index = (day_b - 1)*obj.m + grp_b;
                    
                    aa_val = parent(aa_index);
                    ab_val = parent(ab_index);
                    ba_val = parent(ba_index);
                    bb_val = parent(bb_index);
                    
                    a_total = sum(parent((day_a - 1)*obj.m+1:day_a*obj.m));
                    b_total = sum(parent((day_b - 1)*obj.m+1:day_b*obj.m));
                    
                    max_change = min([ab_val, ba_val, a_total - aa_val, b_total - bb_val]);
                    delta = unifrnd(0, max_change);
                    
                    MutatedChildren(i,aa_index) = aa_val + delta;
                    MutatedChildren(i,bb_index) = bb_val + delta;
                    MutatedChildren(i,ab_index) = ab_val - delta;
                    MutatedChildren(i,ba_index) = ba_val - delta;
                end
            end
            
            out = @mut_dist;
        end
        
        function out = get.CreationFcn(obj)
            
            function Population = f(GenomeLength, FitnessFcn, options)
                Population = obj.getRandomOrderStrategies(400);
            end
            
            out = @f;
        end
    end
    
    %% Comutation intensive methods.
    methods(Access = protected)
        
        function [P, perm] = getRandomOrderStrategies(obj, n)
            disp(['      -> Creating ' num2str(n) ' random order strategies...']);
            P = zeros(n, obj.nvars);
            perm = zeros(n, obj.m);
            
            for i = 1:n
                perm(i,:) = randperm(obj.m);
                P(i,:) = obj.stratmat2input(...
                    obj.VaccinationRestriction.generateOrderStrategyMat(...
                        obj.InitialState, ...
                        perm(i,:), ...
                        obj.days ...
                   ) ...
                );
            end
        end
        
        function P = getRandomFullStrategies(obj, n)
            disp(['      -> Creating ' num2str(n) ' random full strategies...']);
            P = zeros(n, obj.nvars);
            
            for i = 1:n
                P(i, :) = obj.stratmat2input( ...
                   obj.VaccinationRestriction.generateFullStrategyMat(...
                        obj.InitialState, ...
                        obj.days ...
                   ) ...
                );
            end
        end
        
        function [P, perm, score] = getBestOrderStrategies(obj, n)
            disp(['      -> Creating ' num2str(n) ' best order strategies...']);
            
            % Nothing needs to be calculated if n=0;
            if n == 0
                P = zeros(0, obj.nvars);
                perm = zeros(0, obj.m);
                score = zeros(0,1);
                return
            end
            
            % Get all permutations.
            [P, perm, score] = obj.getOrderStrategies();
            
            P = P(1:n, :);
            perm = perm(1:n, :);
            score = score(1:n, :);
        end
        
        function [P] = getInitialPopulationMatrix(obj, varargin)
            disp('    -> Creating the population matrix...');
            
            InitBestOrderStrategies = 0;
            InitRandomOrderStrategies = [];
            InitRandomFullStrategies = 0;
            
            for ii = 1:2:length(varargin)
                switch string(varargin{ii})
                    case "InitBestOrderStrategies"
                        InitBestOrderStrategies = varargin{ii + 1};
                    case "InitRandomOrderStrategies"
                        InitRandomOrderStrategies = varargin{ii + 1};
                    case "InitRandomFullStrategies"
                        InitRandomFullStrategies = varargin{ii + 1};
                end
            end
            
            if isempty(InitRandomOrderStrategies)
                if InitBestOrderStrategies + InitRandomFullStrategies <= 0
                    InitRandomOrderStrategies = 200;
                else
                    InitRandomOrderStrategies = 0;
                end
            end
            
            P = [
                obj.getBestOrderStrategies(InitBestOrderStrategies)
                obj.getRandomOrderStrategies(InitRandomOrderStrategies)
                obj.getRandomFullStrategies(InitRandomFullStrategies)
            ];
            disp(size(P));
        end
        
        function out = getOptions(obj, varargin)
            disp('  -> Creating the optimisation problem options...');
            
            function state = bestInGeneration(options, state, flag)
                [~, Indices] = sort(state.Score);
                
                BestVS = obj.input2strategy(state.Population(Indices(1), :));
                BestVS.plotArea();
                title(['Best in Generation ' num2str(state.Generation)]);
            end
            
            function state = worstInGeneration(options, state, flag)
                [~, Indices] = sort(state.Score);
                
                WorstVS = obj.input2strategy(state.Population(Indices(end), :));
                WorstVS.plotArea();
                title(['Worst in Generation ' num2str(state.Generation)]);
            end
            
            % Default values.
            MaxTime = inf;
            MaxStallTime = inf;
            Display = 'final';
            FunctionTolerance = 1e-4;
            FitnessLimit = 0;
            MaxGenerations = 20000;
            MaxStallGenerations = 50;
            PlotFcn = {@gaplotbestf, @bestInGeneration, @worstInGeneration};
            EliteCount = [];
            
            % Properties
            for ii = 1:2:length(varargin)
                switch string(varargin{ii})
                    case "MaxTime"
                        MaxTime = varargin{ii + 1};
                    case "Display"
                        Display = varargin{ii + 1};
                    case "FitnessLimit"
                        FitnessLimit = varargin{ii + 1};
                    case "FunctionTolerance"
                        FunctionTolerance = varargin{ii + 1};
                    case "MaxGenerations"
                        MaxGenerations = varargin{ii + 1};
                    case "MaxStallGenerations"
                        MaxStallGenerations = varargin{ii + 1};
                    case "MaxStallTime"
                        MaxStallTime = varargin{ii + 1};
                    case "PlotFcn"
                        PlotFcn = varargin{ii + 1};
                    case "EliteCount"
                        EliteCount = varargin{ii + 1};
                end
            end
            
            % Computed default values.
            InitialPopulationMatrix = obj.getInitialPopulationMatrix(varargin{:});
            PopulationSize = height(InitialPopulationMatrix);
            if isempty(EliteCount)
                EliteCount = ceil(0.05 * PopulationSize);
            end
           
            
            % Create options.
            out = optimoptions('ga', ...
                'MaxTime', MaxTime, ...
                'MaxStallTime', MaxStallTime, ...
                'Display', Display, ...
                'UseVectorized', true, ...
                'PopulationSize', PopulationSize, ...
                'ConstraintTolerance', obj.VaccinationRestriction.ConstraintTolerance, ...
                'InitialPopulationMatrix', InitialPopulationMatrix, ...
                'FitnessLimit', FitnessLimit, ...
                'FunctionTolerance', FunctionTolerance, ...
                'MaxGenerations', MaxGenerations, ...
                'CrossoverFcn', {@crossoverintermediate, 1}, ...
                ...'CreationFcn', obj.CreationFcn, ...
                'MutationFcn', obj.MutationFcn, ...
                'MaxStallGenerations', MaxStallGenerations, ...
                'EliteCount', EliteCount, ...
                'PopulationType', 'doubleVector', ...
                'PlotFcn', PlotFcn ...
            );
        
            % Display options.
            disp('[X]___ Optimisation options have been created!___________________________________');
            disp(out);
            disp('_________________________________________________________________________________');
        end
    end
    
    %% Restrictions
    methods
        function [A,b] = solutionSpace(obj)
            A = obj.Aineq;
            b = obj.Bineq;
        end
    end
    
    %% Helper functions
    methods
        
        function M = input2stratmat(obj, X)
            assert(length(X) == obj.nvars);
            M = reshape(X,obj.m,obj.days);
        end
        
        function X = stratmat2input(obj, M)
            X = reshape(M,obj.nvars,1);
        end
        
        function VS = input2strategy(obj, X)
            VS = lib.classes.VaccinationStrategy( ...
                obj.InitialState, ...
                obj.input2stratmat(X) ...
            );
        end
        
        function X = strategy2input(obj, VS)
            X = obj.stratmat2input(VS.Rho);
        end
    end
    
    %% Optimisation function
    methods
        
        function out = getProblem(obj, varargin)
            disp('-> Creating the optimisation problem parameters...');
            
            out = struct( ...
                'fitnessfcn', obj.fitnessfcn, ...
                'nvars', obj.nvars, ...
                'Aineq', obj.Aineq, ...
                'Bineq', obj.Bineq, ...
                'lb', obj.lb, ...
                'solver', 'ga', ...
                'options', obj.getOptions(varargin{:}) ...
            );
            
            disp('[X]___ Optimisation parameters have been created!________________________________');
            disp(out);
            disp('_________________________________________________________________________________');
        end
        
        function [x, perm, score] = getOrderStrategies(obj)
            disp(['        -> Getting all permutations of ' num2str(obj.m) ' groups...']);
            perm = perms(1:obj.m);
            pn = factorial(obj.m);
            % get all order vaccination strategies.
            x = zeros(pn, obj.nvars);
            disp(['        -> Getting all ' num2str(pn) ' order strategies...']);
            for i = 1:pn
                M = obj.VaccinationRestriction.generateOrderStrategyMat(obj.InitialState, perm(i,:), obj.days);
                x(i,:) = obj.stratmat2input(M);
            end
            % Compute and order the scores of those strategies.
            disp(['        -> Compute the scores of the ' num2str(pn) ' order strategies...']);
            unsorted_scores = obj.fitnessfcn(x);
            
            disp(['        -> Order the ' num2str(pn) ' order strategies by their score...']);
            [score, I] = sort(unsorted_scores);
            perm = perm(I,:);
            x = x(I,:);
        end
        
        function [res, fval, x] = runOrderStrategy(obj, perm)
            
            if nargin < 2 || isempty(perm)
                perm = flip(1:obj.m);
            end
            
            res = obj.VaccinationRestriction.generateOrderStrategy(obj.InitialState, perm, obj.days);
            x = obj.strategy2input(res)';
            fval = obj.fitnessfcn(x);
        end
        
        function [res, fval, x] = runRandomStrategy(obj)
            
            res = obj.VaccinationRestriction.generateFullStrategy(obj.InitialState, obj.days);
            x = obj.strategy2input(res)';
            fval = obj.fitnessfcn(x);
        end
        
        function [res, fval, x] = runUniformStrategy(obj)
            
            res = obj.VaccinationRestriction.generateUniformStrategy(obj.InitialState, obj.days);
            x = obj.strategy2input(res)';
            fval = obj.fitnessfcn(x);
        end
        
        function [res, fval, x] = runEmptyStrategy(obj)
            x = zeros(1, obj.nvars);
            fval = obj.fitnessfcn(x);
            res = obj.input2strategy(x);
        end
        
        function [res, fval, perm, x] = runBestOrder(obj)
            disp('=> Getting the best order strategy.');
            
            [x, perm, fval] = obj.getOrderStrategies();
            
            x = x(1,:);
            perm = perm(1,:);
            fval = fval(1);
            
            res = obj.input2strategy(x);
            
        end
        
        function [res, fval, exitflag, output, population, scores, x] = runGeneticAlgorithm(obj, varargin)
            [x, fval, exitflag, output, population, scores] = ga(obj.getProblem(varargin{:}));
            res = lib.classes.VaccinationStrategy( ...
                obj.InitialState, ...
                obj.input2stratmat(x) ...
            );
        end
        
        function [VS, score] = runDifferentialEvolution(obj, varargin)
            [VS, score] = lib.optimize.per_day_optimize(obj, varargin{:});
        end
    end
    
    %% Getting the state.
    methods
        
        function MR = getResult(obj, VS)
            if nargin < 2
                VS = [];
            end
            
            MR = lib.SIRV_model(obj.DeltaT, ...
                'Method', obj.Method, ...
                'Steps', obj.Steps, ...
                'StartDate', obj.StartDate, ...
                'InitialState', obj.InitialState, ...
                'VaccinationStrategy', VS ...
            );
        end
        
    end
    
    %% Plots
    properties(Dependent)
        CostLabel
    end
    
    methods
        function out = get.CostLabel(obj)
            switch string(obj.CostFunc)
                case "total_deaths"
                    out = "Total Amount of Deaths";
                otherwise
                    out = "";
            end
        end
    end
end

