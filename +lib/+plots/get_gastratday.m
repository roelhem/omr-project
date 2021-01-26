function [out, day_indices] = get_gastratday(m, day, th_max)
%GET_GASTRATDAY Summary of this function goes here
%   Detailed explanation goes here
    
    if nargin < 1
        m = [];
    end
    
    if nargin < 2
        day = 1;
    end
    
    day_indices = [];
    if ~isempty(m)
        day_indices = (day - 1)*m+1:day*m;
    end

    function state = gastratday(options, state, flag)
        [~, I] = sort(state.Score);
        
        best = state.Population(I(1), day_indices);
        lib.plots.stratday(best, th_max);
        title('Best in generation ');
    end

    out = @gastratday;

end

