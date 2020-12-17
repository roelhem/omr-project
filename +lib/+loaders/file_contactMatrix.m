function [T, data, groups, groupBoundaries] = file_contactMatrix(filename)
%FILE_CONTACTMATRIX Reads and validates the date from a contact matrix csv-file.


%% Get the right file path
filepath = filename;
if ~isfile(filepath)
    filepath = strcat('data/contact_matrices/', filename);
end
assert(isfile(filepath), ['Could not find file at path "' filepath '".'])

%% Loading the data from the csv-file and assert that the format is correct.
data = readmatrix(filepath, 'Range', [2, 1]);

n = width(data);
assert(height(data) == n, ['The matrix from "' filepath '" is not a square matrix.']);

%% Loading the groups from the matrix.
groups = readcell(filepath, 'Range', [1 1 1 n], 'DatetimeType', 'text')';
groupBoundaries = lib.utils.strToBoundaries(groups);
groups = lib.utils.boundariesToCat(groupBoundaries);

%% Put everything in a table.
T = array2table(data);
T.Properties.VariableNames = string(groups);
T.Properties.RowNames = string(groups);

Group = table(groups, groupBoundaries, 'VariableNames', {'Category', 'Boundaries'});

T = addvars(T, Group);


end