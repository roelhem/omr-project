function [S,I,R,V] = SIRV_EulerForward(Szero,Izero,Rzero,Vzero,beta,alpha,nu,dt,n)
%SIR_EULERFORWARD Runs a simple SIR-model using the Euler Forward method.
%   It will always start at `t=0`.

%% Get the amount of age groups and total amount of people.
m = length(Szero);
N = (Szero + Izero + Rzero+Vzero)';

%% Ensure that the input has the right demensions.
% Checking the initial values.
assert(length(Izero) == m, 'Izero has incompatible length')
assert(length(Rzero) == m, 'Rzero has incompatible amount of groups.')
assert(length(Vzero) == m, 'Vzero has incompatible amount of groups.')
% Converting the parameters.
beta  = lib.utils.asParamMat(beta , m);
alpha = lib.utils.asParamMat(alpha, m);
nu    = lib.utils.asParamMat(nu   , m);

%% Allocating memory for the result.
S = zeros(m, n);
I = zeros(m, n);
R = zeros(m, n);
V = zeros(m, n);

%% Initial conditions.
S(:,1) = Szero;
I(:,1) = Izero;
R(:,1) = Rzero;
V(:,1) = Vzero;


%% Reccursively computing the values.
for i = 1:n-1
    % Getting the transmission between groups.
    StoI = (beta * (I(:,i) ./ N)) .* S(:,i);
    ItoR = alpha * I(:,i);
    StoV = nu * S(:,i);

    % Reccursively update the state.
    S(:,i+1) = S(:,i) + (- StoI - StoV) * dt;
    I(:,i+1) = I(:,i) + (  StoI - ItoR) * dt;
    R(:,i+1) = R(:,i) + (         ItoR) * dt;    
    V(:,i+1) = V(:,i) + (         StoV) * dt;
end
end

