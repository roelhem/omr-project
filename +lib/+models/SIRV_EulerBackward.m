function [S,I,R,V] = SIRV_EulerBackward(Szero,Izero,Rzero,Vzero,beta,alpha,nu,dt,n)
%SIRV_EULERFORWARD Runs the SIRV-model using the Euler Backward method.
%   It will always start at `t=0`.

%% Get the amount of age groups and total amount of people.
m = height(Szero);
N = (Szero + Izero + Rzero + Vzero);

%% Ensure that the input has the right demensions.
% Checking the initial values.
assert(length(Izero) == m, 'Izero has incompatible length');
assert(length(Rzero) == m, 'Rzero has incompatible amount of groups.');
assert(length(Vzero) == m, 'Vzero has incompatible amount of groups.');
% Converting the parameters.
beta  = lib.utils.asParamMat(beta, m);
alpha = lib.utils.asParamMat(alpha, m);
% Checking the nu function.
assert(height(nu) == m, 'nu has incompatible amount of groups.');
assert(width(nu) == n, 'nu has incompatible amount of steps.');

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
    StoI = @(x_I, x_S) (beta * x_I ./ N) .* x_S;
    ItoR = @(x_I) alpha * x_I;
    StoV = @(i, x_S) (nu(:, i) .* x_S);
    ItoV = @(i, x_I) (nu(:, i) .* x_I);
    RtoV = @(i, x_R) (nu(:, i) .* x_R);

    % Get the optimistic values.
    S_new = S(:,i) + (- StoI(I(:,i), S(:,i))                - StoV(i, S(:,i)) ) * dt;
    I_new = I(:,i) + (  StoI(I(:,i), S(:,i)) - ItoR(I(:,i)) - ItoV(i, I(:,i)) ) * dt;
    R_new = R(:,i) + (                         ItoR(I(:,i)) - RtoV(i, R(:,i)) ) * dt;
    
    % Set the new values.
    S(:,i+1) = S(:,i) + (- StoI(I_new, S_new)               - StoV(i+1, S_new) ) * dt;
    I(:,i+1) = I(:,i) + (  StoI(I_new, S_new) - ItoR(I_new) - ItoV(i+1, I_new) ) * dt;
    R(:,i+1) = R(:,i) + (                       ItoR(I_new) - RtoV(i+1, R_new) ) * dt;
end

end

