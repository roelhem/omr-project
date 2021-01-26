function out = get_Reff_vectorized(InitialState)
%GET_TOTALDEATHS_VECTORIZED Calculates the total deaths of multiple
% vaccination plans.

    %% AHEAD OF TIME constant computation.
    % Array sizes and indices.
    m = InitialState.m;
    
    % Population vector.
    N = InitialState.N;
    N_tot = sum(N);
    ReprNum = InitialState.ReprNum;
    
    % Matrix to determine StoI
    B = (InitialState.Beta * ((1 ./ N) .* eye(m)))';
    % Matrix to determine ItoR
    A = InitialState.Alpha';
    
    % Initial values;
    S_initial = InitialState.S';
    I_initial = InitialState.I';
    U_initial = InitialState.U';
    
    %% AHEAD OF TIME parameter validation.
    
    %% AT TIME COMPUTATIONS.
    % Defining the vectorized cost function. This should run as fast as
    % possible.
    function cost = Reff_vectorized(x)
        % Get the genetic algorithm population size.
        p = height(x);
        
        % Initialize the cost.
        cost = zeros(p, 1);
        
        % Initialize initial values.
        S = repmat(S_initial, p, 1);
        I = repmat(I_initial, p, 1);
        U = repmat(U_initial, p, 1);
        
        % Compute cost using EulerForward.
        for li = 1:m:width(x)
            % Reusable variables.
            Rho = x(:, li:li+m-1);
            Nu = Rho ./ U;
            StoI = (I * B) .* S;
            
            % Increment the cost.
            cost = cost + ReprNum * sum(S, 2) / N_tot;
            
            % Computing S, I and U for the next step.
            S = S - StoI            - Nu .* S;
            I = I + StoI - (A .* I) - Nu .* I;
            U = U - Rho;
        end
    end

    %% RESULT
    % Returning the function handle.
    out = @Reff_vectorized;

end

