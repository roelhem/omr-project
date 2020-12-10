function result = fixDateAutocorrect(input)
%FIXDATEAUTOCORRECT Replaces date-names with date-numbers in a string.

% Defining which values should map to which numbers.
DateMap = containers.Map([
    {'Jan'}
    {'Feb'}
    {'Mar'}
    {'Apr'}
    {'May'}
    {'Jun'}
    {'Jul'}
    {'Aug'}
    {'Sep'}
    {'Oct'}
    {'Nov'}
    {'Dec'}
], (1:12)');

% Initialising the result.
result = input;

% Replace each month name with the value from the DateMap.
for month = keys(DateMap)
    result = strrep(result, month, num2str(DateMap(month{1})));
end


end

