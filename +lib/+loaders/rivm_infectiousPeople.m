function result = rivm_infectiousPeople()
%RIVM_ Summary of this function goes here
%   Source: https://data.rivm.nl/covid-19/COVID-19_prevalentie.json

% Get the data-response from the website of RIVM.
response = webread('https://data.rivm.nl/covid-19/COVID-19_prevalentie.json');

n = height(response);

result = table('Size', [n 5],...
               'VariableTypes', {'datetime';'doubleNaN';'doubleNaN';'doubleNaN';'categorical'},...
               'VariableNames', {'Date';'prev_low';'prev_avg';'prev_up';'population'}...
         );

for i = 1:n
    result.Date(i) = response{i,1}.Date;
    if isfield(response{i,1}, 'prev_low')
        result.prev_low(i) = response{i,1}.prev_low;
    end
    
    if isfield(response{i,1}, 'prev_up')
        result.prev_up(i)  = response{i,1}.prev_up;
    end
    
    if isfield(response{i,1}, 'prev_avg')
        result.prev_avg(i)  = response{i,1}.prev_avg;
    elseif isfield(response{i,1}, 'prev_low') && isfield(response{i,1}, 'prev_up')
        result.prev_avg(i)  = (response{i,1}.prev_low + response{i,1}.prev_up) / 2;
    end
    
    result.population(i) = response{i,1}.population;
end

end

