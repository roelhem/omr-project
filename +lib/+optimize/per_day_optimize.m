function [VS, score] = per_day_optimize(OptProblem, varargin)
%PER_DAY_OPTIMIZE Optimise the problem per day.

    %% Get ahead of time variables.
    m = OptProblem.m;
    base_fitnessfcn = OptProblem.fitnessfcn;
    
    %% Get tweakable parameters.
    Iterations = 1000;
    PopulationSize = 50;
    CrossoverProbability = 0.2;
    ConstraintTolerance = OptProblem.VaccinationRestriction.ConstraintTolerance;
    
    for ii = 1:2:length(varargin)
        switch string(varargin{ii})
            case "Iterations" % Amount of iterations.
                Iterations = varargin{ii + 1};
            case "PopulationSize" % Size of the algorithm population.
                PopulationSize = varargin{ii + 1};
            case "CrossoverProbability" % Probability that crossover will be applied.
                CrossoverProbability = varargin{ii + 1};
            case "ConstraintTolerance"
                ConstraintTolerance = varargin{ii + 1};
        end
    end
    
    %% Initialisation.
    
    disp('=> INITIALISING THE DIFFERENTIAL EVOLUTION ALGORITM...');
    
    % Aliases
    n = OptProblem.days;
    p = PopulationSize;
    tol = ConstraintTolerance;
    
    % Population
    Population = zeros(p, m);
    Scores = inf(p, 1);
    
    % Changing variables.
    result_input = zeros(1, m*n);
    score = inf;
    N_left = OptProblem.N * OptProblem.VaccinationRestriction.Efficacy;
    N_left = N_left';
    
    %% Plot Setup
    figure;
    
    subplot(2,1,1);
    bar(Population, 'stacked');
    ylim([0 OptProblem.VaccinationRestriction.EffMaxPerDay]);
    xlabel('Best to worst strategy');
    ylabel('Vaccinations \rho');
    
    subplot(2,1,2);
    area(zeros(n,m));
    ylim([0 OptProblem.VaccinationRestriction.EffMaxPerDay]);
    xlim([0 n]);
    ylabel('Vaccinations \rho');
    xlabel('Days');
    
    drawnow;
    
    
    %% Run the differential evolution for each day.
    for day = 1:n
        
        % Output the progress.
        disp(['=> OPTIMISING DAY ' num2str(day)]);
        
        % Compute the amount of available vaccines.
        maxvacc = min(OptProblem.VaccinationRestriction.EffMaxPerDay, sum(N_left));
        if maxvacc < tol
            disp([' --> [TERMINADED] All available vaccines used. ' num2str(day)]);
            break;
        end
        
        %% Differential evolution algorithm
        
        % Get the fitness function.
        [fitnessfcn, change_indices] = lib.utils.getPartialFitnessFunc([result_input, zeros(1,m*40)], day, m, base_fitnessfcn);
        
        % Get the initial population.
        Population(:,:) = lib.utils.randomVaccDist(N_left, maxvacc, p, tol);
        Scores(:) = fitnessfcn(Population);
        
        % Run each iteration.
        for i = 1:Iterations
            
            % MUTATION
            beta = unifrnd(0.2, 0.8, p, m);
            picked = randPermMat(p, 3);
            y = Population(picked(:,1)) + beta .* (Population(picked(:,2)) - Population(picked(:,3)));
            y = max(0, y);
            y = min(N_left, y);
            
            % CROSSOVER
            z = Population;
            j = repmat(1:m, p, 1) == randi([1 m], p, 1);
            z(j) = y(j);
            z = z .* (maxvacc ./ max(0, sum(z, 2)));
            
            % SELECTION
            z_scores = fitnessfcn(z);
            better = z_scores < Scores;
            Scores(better) = z_scores(better);
            Population(better, :) = z(better, :);
            
            % Plotting
             if mod(i, 1000) == 0
                 subplot(2,1,1);
                 [~, I] = sort(Scores, 'descend');
                 bar(Population(I,:), 'stacked');
                 title(['Strategies for day ' num2str(day) ' at iteration ' num2str(i)]);
                 
                 drawnow;
             end
        end
        
        %% Get intermediate result.
        [score, I] = min(Scores);
        best = Population(I, :);
        result_input(change_indices) = best;
        Rho = OptProblem.input2stratmat(result_input);
        
        if isnan(score)
            break;
        end
        
        %% Determine next values
        N_left = max(0, N_left - best);
        
        %% Plot results.
        subplot(2,1,2);
        area(Rho');
        title(['Result till day ' num2str(day) ', cost ' num2str(score)]);
        
        drawnow;
        
    end
    
    %% Get the result.
    VS = OptProblem.input2strategy(result_input);

end

function out = randPermMat(p, amount)
    out = zeros(p, amount);
    for pi = 1:p
        out(pi,:) = randperm(p, amount);
    end
end
