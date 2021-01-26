function out = get_maxI_vectorized(InitialState)
%GET_MAXI_VECTORIZED Calculates the total deaths of multiple
% vaccination plans.

    %% AHEAD OF TIME constant computation.
    % Array sizes and indices.
    m = InitialState.m;
    
    % Population vector.
    N = InitialState.N;
    
    % Matrix to determine StoI
    B = (InitialState.Beta * ((1 ./ N) .* eye(m)))';
    % Matrix to determine ItoR
    A = InitialState.Alpha';
    
    % Initial values;
    S_initial = InitialState.S';
    I_initial = InitialState.I';
    U_initial = InitialState.U';
    
    %% AT TIME COMPUTATIONS.
    % Defining the vectorized cost function. This should run as fast as
    % possible.
    function cost = maxI_vectorized(x)
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
            I_sum = sum(I, 2);
            Larger = I_sum > cost;
            cost(Larger) = I_sum(Larger);
            
            % Computing S, I and U for the next step.
            S = S - StoI            - Nu .* S;
            I = I + StoI - (A .* I) - Nu .* I;
            U = U - Rho;
        end
    end

    %% RESULT
    % Returning the function handle.
    out = @maxI_vectorized;

end

