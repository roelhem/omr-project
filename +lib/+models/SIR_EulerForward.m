function [S,I,R] = SIR_EulerForward(Szero,Izero,Rzero,beta,alpha,dt,n)
%SIR_EULERFORWARD Runs a simple SIR-model using the Euler Forward method.
%   It will always start at `t=0`.

%% Get the amount of age groups and total amount of people.
disp(Szero);
m = height(Szero);
disp(m);
N = (Szero + Izero + Rzero);

%% Ensure that the input has the right demensions.
% Checking the initial values.
assert(length(Izero) == m, 'Izero has incompatible length')
assert(length(Rzero) == m, 'Rzero has incompatible amount of groups.')
% Converting the parameters.
beta  = lib.utils.asParamMat(beta, m);
alpha = lib.utils.asParamMat(alpha, m);

%% Allocating memory for the result.
S = zeros(m, n);
I = zeros(m, n);
R = zeros(m, n);

%% Initial conditions.
S(:,1) = Szero;
I(:,1) = Izero;
R(:,1) = Rzero;


%% Reccursively computing the values.
for i = 1:n-1
    % Getting the transmission between groups.
    StoI = (beta * (I(:,i) ./ N)) .* S(:,i);
    ItoR = alpha * I(:,i);

    % Reccursively update the state.
    S(:,i+1) = S(:,i) + (- StoI       ) * dt;
    I(:,i+1) = I(:,i) + (  StoI - ItoR) * dt;
    R(:,i+1) = R(:,i) + (         ItoR) * dt;
end
end

