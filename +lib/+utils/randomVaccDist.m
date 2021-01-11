function [out, dist] = randomVaccDist(N, maxvacc, p, tol)
%RANDOMVACCDIST Summary of this function goes here
%   Detailed explanation goes here

%% Get parameters.
% Tolerance.
if nargin < 4
    tol = 1e-5;
end
tol = abs(tol);

% Get default max value.
if nargin < 2 || isempty(maxvacc)
    maxvacc = sum(N);
else
    maxvacc = min(maxvacc, sum(N));
end

% Get default population size.
if nargin < 3 || isempty(p)
    p = 1;
end

% Format N
if width(N) == 1
    N = N';
end
assert(height(N) == 1, 'N has to have width (or height) of size 1');

% Get the amount of age groups.
m = width(N);

%% Initialise result
out = zeros(p, m);

%% Handle edge cases.
% Return zeros when max <= 0. (1e-5 is the tolerance)
if maxvacc < tol
     return
end

%% Create the output.
underflow    = repmat(N, p, 1);
target_total = repmat(maxvacc, p, 1);
corr         = ones(p, 1, 'logical');
while any(corr) && any(underflow >= tol, 'all')
    % Get the target total
    target_total = min(sum(underflow, 2), target_total);
    
    % Distribute the population randomly.
    add = rand(height(underflow), m) .* underflow;
    add = add .* (target_total ./ sum(add, 2));
    
    % Correct the overflow by setting to the maximum value.
    pop_overflow = add > underflow;
    add(pop_overflow) = underflow(pop_overflow);
    
    % Update the result.
    out(corr,:) = out(corr,:) + add;
    
    % Determine the values that haven't meet there target total.
    add_total = sum(add, 2);
    target_total = target_total - add_total;
    review = target_total >= tol;
    
    % Determine values for next iteration.
    underflow = underflow(review, :) - add(review, :);
    target_total = target_total(review);
    corr(corr) = review;
end


end

