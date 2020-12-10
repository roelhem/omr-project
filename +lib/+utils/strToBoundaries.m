function result = strToBoundaries(input)
% STRTOBOUNDARIES Extracts the boundaries from a string.

    % Ensure that the input is a string.
    input = lib.utils.fixDateAutocorrect(string(input));
    
    % Assert that the width of the input is 1.
    assert(width(input) == 1, "Input width is not equal to 1.")
    
    % Function to parse one string value.
    function boundaries = parseStr(val)
        if(contains(val, "-"))
            parts = str2double(split(val, "-"));
            boundaries = [min(parts), max(parts)];
        elseif(contains(val, "+"))
            lowerBound = str2double(extractBefore(val, "+"));
            boundaries = [lowerBound, inf];
        elseif(contains(val, "<"))
            upperBound = str2double(extractAfter(val, "<"));
            boundaries = [0, upperBound];
        else
            error(strcat("The string '", val, "' cannot be converted to a boundary."));
        end
    end

    % Parse each string in the array.
    result = cell2mat(arrayfun(@parseStr, input, 'UniformOutput', false));
end