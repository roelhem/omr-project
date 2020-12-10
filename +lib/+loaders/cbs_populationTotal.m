function output = cbs_populationTotal()
%CBS_POPULATIONTOTAL Gets the total population from the CBS website.
%   Detailed explanation goes here

    %% Get the url
    endpoint = "https://opendata.cbs.nl/ODataApi/odata/84646NED/TypedDataSet";
    filter   = "substring(Perioden,0,4) eq '2020' and substring(Perioden,4,2) eq 'JJ'";
    select   = "Bevolking_1,Leeftijd,Geslacht";

    url = strcat(endpoint, "?$filter=", filter, "&$select=", select);

    %% Request the meaning of the values in the dataset.
    % Age values.
    age_response = webread('https://opendata.cbs.nl/ODataApi/odata/84646NED/Leeftijd?$select=Key,Title,CategoryGroupID');
    age_table = struct2table(age_response.value);
    age2title_map = containers.Map(age_table.Key, age_table.Title);
    age2category_map = containers.Map(age_table.Key, age_table.CategoryGroupID);

    % Gender values.
    %gender_response = webread('https://opendata.cbs.nl/ODataApi/odata/84646NED/Geslacht?$select=Key,Title');
    %gender_table = struct2table(gender_response.value);
    gender_map = containers.Map([
        "3000   "
        "4000   "
        "T001038"
    ],[
        {'Male'}
        {'Female'}
        {'Total'}
    ]);

    %% Request the data.
    response = webread(url);
    T = struct2table(response.value);

    %% Format the data.
    T.Leeftijd = categorical(T.Leeftijd);
    %output = arrayfun(@(x)age2category_map(string(x)), table.Leeftijd, 'UniformOutput', false);
    T.Geslacht = categorical(arrayfun(@(x)gender_map(string(x)), T.Geslacht, 'UniformOutput', false));
    T = unstack(T, 'Bevolking_1', 'Geslacht');

    % Add the names of the age groups.
    T.AgeGroupTitle = arrayfun(@(x)age2title_map(string(x)), T.Leeftijd, 'UniformOutput', false);

    % Add the categories (type of age group) of the age groups.
    categoryIds = [1,(2:11),13,14,(15:16)];
    categoryNames = repelem(["Total","One","Five","Ten","Legal"],[1,10,1,1,2]);
    T.GroupCategoryId = arrayfun(@(x)age2category_map(string(x)), T.Leeftijd);
    T.GroupCategory = categorical(T.GroupCategoryId, categoryIds, categoryNames);
    
    % Get the boundaries of the group from the group titles.
    regex = "^Totaal leeftijd$|^(?<Exactly>\d+) jaar$|^(?<LowerBound>\d+) tot (?<UpperBound>\d+) jaar$|^(?<LowerBound>\d+) jaar of ouder$";
    matches = struct2table(cell2mat(regexp(string(T.AgeGroupTitle), regex, 'names')));
    
    function [lowerBound, upperBound] = getBound(exact, low, high)
        if(~isempty(exact{1}))
            lowerBound = str2double(exact{1});
            upperBound = str2double(exact{1});
        elseif(~isempty(high{1}) && ~isempty(low{1}))
            lowerBound = str2double(low{1});
            upperBound = str2double(high{1});
        elseif(~isempty(low{1}))
            lowerBound = str2double(low{1});
            upperBound = inf;
        else
            lowerBound = 0;
            upperBound = inf;
        end
    end

    T.GroupBoundaries = table2array(rowfun(@getBound, matches, 'OutputVariableNames', {'LowerBound' 'UpperBound'}));
    
    Group = T(:, {'GroupBoundaries'});
    Group.Properties.VariableNames = {'Boundaries'};
    Group = addvars(Group, lib.utils.boundariesToCat(Group.Boundaries), 'Before', 'Boundaries','NewVariableName', 'Category');
    
    TResultTable = T(:, {'Male','Female','Total'});
    TResultTable = addvars(TResultTable, Group, 'Before', 'Male');
    TResultTable.Properties.RowNames = T.AgeGroupTitle;
    
    %% Set Data Display.
     output.Total = TResultTable(T.GroupCategory == "Total",:);
     output.One = TResultTable(T.GroupCategory == "One",:);
     output.Five = TResultTable(T.GroupCategory == "Five",:);
     output.Ten = TResultTable(T.GroupCategory == "Ten",:);
     output.Legal = TResultTable(T.GroupCategory == "Legal",:);

end

