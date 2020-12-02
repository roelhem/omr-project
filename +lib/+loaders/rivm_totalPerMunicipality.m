function result = rivm_totalPerMunicipality()
%GET_COVIDNUMBERPERMUNICIPALITY Gets the data per municipality. 
%   Source: https://data.rivm.nl/covid-19/COVID-19_aantallen_gemeente_cumulatief.json

% Get the data-response from the website of RIVM.
response = webread('https://data.rivm.nl/covid-19/COVID-19_aantallen_gemeente_cumulatief.json');

% Convert the data to a table.
result = struct2table(response);

% Format the results.
result.Date_of_report = datetime(result.Date_of_report);
result.Municipality_code = categorical(lib.utils.ensureChar(result.Municipality_code, {'Unknown'}));

notChars = cellfun(@(x)(not(ischar(x))), result.Municipality_name);
result.Municipality_name(notChars) = strcat('Unknown[', result.Province(notChars),']');
result.Municipality_name = categorical(result.Municipality_name);
result.Province = categorical(result.Province);

end

