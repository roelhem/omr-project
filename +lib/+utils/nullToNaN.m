function cells = nullToNaN(cells)
%NULLTONAN Summary of this function goes here
%   Detailed explanation goes here

    function y = convert(x)
        if isempty(x)
            y = NaN;
        else
            y = x(1);
        end
    end
    cells = cellfun(@convert, cells);
end

