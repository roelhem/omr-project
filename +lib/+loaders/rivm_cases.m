function result = rivm_cases()
%GET_COVIDNUMBERPERMUNICIPALITY Gets all reported covid cases in The Netherlands. 
%   Source: https://data.rivm.nl/covid-19/COVID-19_casus_landelijk.json

% Get the data-response from the website of RIVM.
response = webread('https://data.rivm.nl/covid-19/COVID-19_casus_landelijk.json');

% Convert the data to a table.
result = struct2table(response);
result.Date_file            = datetime(result.Date_file);
result.Date_statistics      = datetime(result.Date_statistics);
result.Date_statistics_type = categorical(result.Date_statistics_type);
result.Agegroup             = categorical(result.Agegroup);
result.Sex                  = categorical(result.Sex);
result.Hospital_admission(cellfun(@(x)(not(ischar(x))), result.Hospital_admission)) = {'Unknown'};
result.Hospital_admission   = categorical(result.Hospital_admission);
result.Deceased(cellfun(@(x)(not(ischar(x))), result.Deceased)) = {'Unknown'};
result.Deceased             = categorical(result.Deceased);
result.Province(cellfun(@(x)(not(ischar(x))), result.Province)) = {'Unknown'};
result.Province = categorical(result.Province);
result.Municipal_health_service(cellfun(@(x)(not(ischar(x))), result.Municipal_health_service)) = {'Unknown'};
result.Municipal_health_service = categorical(result.Municipal_health_service);

I = not(cellfun(@isempty, result.Week_of_death));
weekStrings = cell2mat(result.Week_of_death(I));
years = str2num(weekStrings(:,1:4));
weeks = str2num(weekStrings(:,5:6));

Week_of_death = NaN(height(result),1);
Week_of_death(I) = weeks;
result.Week_of_death = Week_of_death;

Year_of_death = NaN(height(result),1);
Year_of_death(I) = years;
result = addvars(result, Year_of_death, 'Before', 'Week_of_death');

end
