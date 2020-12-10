function output = boundariesToCat(input)
%BOUNDARIESTOSTR Takes a matrix of boundaries and returns the categories to represent them.

    output = strings(height(input), 1);
    
    for i = 1:height(input)
        if isinf(input(i,2))
            output(i) = strcat(num2str(input(i,1)), "+");
        else
            output(i) = strcat(num2str(input(i,1)), "-", num2str(input(i,2)));
        end
    end

    output = categorical(output);
end

