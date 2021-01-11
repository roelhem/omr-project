function [out, change_indices, x_template] = getPartialFitnessFunc(x_template, day, m, fitnessfcn)
%PARTIALFITNESSFUNC Returns a fitness function that only takes the input of
%one day.

    change_indices = (day - 1)*m+1:day*m;
    
    if width(x_template) == 1
        x_template = x_template';
    end

    function costs = partial_fitnessfcn(x)
        new_x = repmat(x_template, height(x), 1);
        new_x(:, change_indices) = x;
        costs = fitnessfcn(new_x);
    end

    out = @partial_fitnessfcn;

end

