function out = SIRV_EulerForward_vectorized(InitialState)
%SIRV_EULERFORWARD_VECTORIZED Summary of this function goes here
%   Detailed explanation goes here

%% AHEAD OF TIME constant computation.
    % Array sizes and indices.
    m = InitialState.m;
    
    % Population vector.
    N = InitialState.N;
    
    % Matrix to determine StoI
    B = (InitialState.Beta * ((1 ./ N) .* eye(m)))';
    % Matrix to determine ItoR
    A = InitialState.Alpha(1);
    
    % Initial values;
    S_initial = InitialState.S';
    I_initial = InitialState.I';
    U_initial = InitialState.U';
    
    %% AT TIME COMPUTATIONS.
    % Defining the vectorized cost function. This should run as fast as
    % possible.
    function [S,I,U,StoI,Nu] = f(x)
        
        % Get the genetic algorithm population size.
        p = height(x);
        l = length(x);
        n = l / m;
        
        % Initialize the values.
        S = zeros(p, m, n);
        I = zeros(p, m, n);
        U = zeros(p, m, n);
        StoI = zeros(p, m, n);
        Nu = zeros(p, m, n);
        
        % Initialize initial values.
        S(:,:,1) = repmat(S_initial, p, 1);
        I(:,:,1) = repmat(I_initial, p, 1);
        U(:,:,1) = repmat(U_initial, p, 1);
        
        % Compute cost using EulerForward.
        for i = 1:n-1
            % Reusable variables.
            Rho = x(:, (i - 1)*m + 1:i*m);
            Nu(:,:,i) = Rho ./ U(:,:,i);
            % Computing StoI, S, I and U for the next step.
            StoI(:,:,i) = (I(:,:,i) * B) .* S(:,:,i);
            S(:,:,i+1) = S(:,:,i) - StoI(:,:,i)                   - Nu(:,:,i) .* S(:,:,i);
            I(:,:,i+1) = I(:,:,i) + StoI(:,:,i) - (I(:,:,i) * A)  - Nu(:,:,i) .* I(:,:,i);
            U(:,:,i+1) = U(:,:,i) - Rho;
        end
        
        % Compute last StoI.
        StoI(:,:,n) = (I(:,:,n) * B) .* S(:,:,n);
    end

    %% Returning the result
    out = @f;

end

