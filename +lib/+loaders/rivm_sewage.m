function result = rivm_sewage()
%RI Summary of this function goes here
%   Source: https://data.rivm.nl/covid-19/COVID-19_rioolwaterdata.json

% Get the data-response from the website of RIVM.
response = webread('https://data.rivm.nl/covid-19/COVID-19_rioolwaterdata.json');

% Convert the data to a table.
result = struct2table(response);


result.Date_measurement = datetime(result.Date_measurement);
result.RWZI_AWZI_name = categorical(result.RWZI_AWZI_name);
result.Security_region_code = categorical(result.Security_region_code);
result.Security_region_name = categorical(result.Security_region_name);
result.Percentage_in_security_region = str2double(strrep(result.Percentage_in_security_region, ',', '.'));
result.RNA_per_ml = lib.utils.nullToNaN(result.RNA_per_ml);
result.RNA_flow_per_100000 = lib.utils.nullToNaN(result.RNA_flow_per_100000);

end

