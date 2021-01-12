function [res, fval, exitflag, output, population, scores, x] = one_group_per_day_optimize(OptProblem,varargin)
%ONE_GROUP_PER_DAY_OPTIMIZE Summary of this function goes here
%   Detailed explanation goes here

    %% Get parameters.
    PopulationSize = 50;
    
    for ii = 1:2:length(varargin)
        switch string(varargin{ii})
            case "PopulationSize" % Size of the algorithm population.
                PopulationSize = varargin{ii + 1};
        end
    end


    %% Ahead of time computations.
    N_days = floor(OptProblem.N / OptProblem.VaccinationRestriction.MaxPerDay);
    
    maxvacc = OptProblem.VaccinationRestriction.EffMaxPerDay;
    n = OptProblem.days;
    m = OptProblem.m;
    p = PopulationSize;

    %% Get the I_map and Generate random permutations.
    I_map = zeros(1, n);
    l = 0;
    for mi = 1:m
        a = N_days(mi);
        I_map(l + (1:a)) = mi;
        l = l + a;
    end
    nvars = l*m;

    InitialPopulation = zeros(p, nvars);
    for pi = 1:p
        rand_perm = I_map(randperm(l));
        for li = 1:l
            InitialPopulation(pi, (li-1)*m + rand_perm(li)) = true;
        end
    end
    
    %% Boundaries
    A = repmat(eye(m),1,l);
    b = N_days;
    
    %% Helper functions
    function out = toInput(x)
        out = zeros(height(x), n*m);
        out(:, 1:width(x)) = x * maxvacc;
    end

    function costs = fitnessfcn(x)
        costs = OptProblem.fitnessfcn(toInput(x));
    end

    function MutatedChildren = MutationFcn(parents, options, nvars, FitnessFcn, state, thisScore, thisPopulation)
        MutatedChildren = thisPopulation(parents,:);
        
        ps = height(parents);
        for psi = 1:ps
            swap_vals = randperm(l,2);
            I_a = 1:m + (swap_vals(1) - 1)*m;
            I_b = 1:m + (swap_vals(2) - 1)*m;
            temp = MutatedChildren(psi, I_a);
            MutatedChildren(psi, I_a) = MutatedChildren(psi, I_b);
            MutatedChildren(psi, I_b) = temp;
        end
        
    end

    
    %% Get genetic algorithm options and properties.
    function state = bestInGeneration(options, state, flag)
        [~, Indices] = sort(state.Score);

        BestVS = OptProblem.input2strategy(toInput(state.Population(Indices(1), :)));
        BestVS.plotArea();
        title(['Best in Generation ' num2str(state.Generation)]);
    end

    function state = worstInGeneration(options, state, flag)
        [~, Indices] = sort(state.Score);

        WorstVS = OptProblem.input2strategy(toInput(state.Population(Indices(end), :)));
        WorstVS.plotArea();
        title(['Worst in Generation ' num2str(state.Generation)]);
    end
    
    options = optimoptions('ga', ...
        'UseVectorized', true, ...
        'InitialPopulationMatrix', InitialPopulation, ...
        'PlotFcn', {@gaplotbestf, @bestInGeneration, @worstInGeneration}, ...
        'MutationFcn', @MutationFcn, ...
        'PopulationType', 'doubleVector', ...
        'PopulationSize', PopulationSize ...
    );

    problem = struct( ...
        'fitnessfcn', @fitnessfcn, ...
        'nvars', nvars, ...
        'solver', 'ga', ...
        'lb', zeros(1,nvars), ...
        'ub', ones(1,nvars), ...
        'Aineq', A, ...
        'Bineq', b, ...
        'IntCon', 1:nvars, ...
        'options', options ...
    );

    %% Run the algoritm 
    [x, fval, exitflag, output, population, scores] = ga(problem);
    res = lib.classes.VaccinationStrategy( ...
        OptProblem.InitialState, ...
        OptProblem.input2stratmat(toInput(x)) ...
    );

end


