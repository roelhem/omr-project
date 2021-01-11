function result = rivm_reproduction()
%RIVM_REPRODUCTION Summary of this function goes here
%   Source: https://data.rivm.nl/covid-19/COVID-19_reproductiegetal.json

% Reproduction Number
response = webread('https://data.rivm.nl/covid-19/COVID-19_reproductiegetal.json');
n = height(response);

result = table('Size', [n 5],...
               'VariableTypes', {'datetime';'doubleNaN';'doubleNaN';'doubleNaN';'categorical'},...
               'VariableNames', {'Date';'Rt_low';'Rt_avg';'Rt_up';'population'}...
         );

for i = 1:n
    result.Date(i) = response{i,1}.Date;
    if isfield(response{i,1}, 'Rt_low')
        result.Rt_low(i) = str2double(response{i,1}.Rt_low);
    end
    
    if isfield(response{i,1}, 'Rt_up')
        result.Rt_up(i)  = str2double(response{i,1}.Rt_up);
    end
    
    if isfield(response{i,1}, 'Rt_avg')
        result.Rt_avg(i)  = str2double(response{i,1}.Rt_avg);
    elseif isfield(response{i,1}, 'Rt_low') && isfield(response{i,1}, 'Rt_up')
        result.Rt_avg(i)  = (str2double(response{i,1}.Rt_up) + str2double(response{i,1}.Rt_low)) / 2;
    end
    
    result.population(i) = response{i,1}.population;
end

end

