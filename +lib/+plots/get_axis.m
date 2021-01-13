function [ax] = get_axis(varargin)
%GET_AXIS Summary of this function goes here
%   Detailed explanation goes here

    Axis = [];
    Figure = [];

    for ii = 1:2:length(varargin)
        switch string(varargin{ii})
            case "Axis"
                Axis = varargin{ii + 1};
            case "Figure"
                Figure = varargin{ii + 1};
        end
    end
    
    if isempty(Axis)
        if isa(Figure, 'matlab.ui.Figure') && ~isempty(Figure.CurrentAxes)
            Axis = Figure.CurrentAxes;
        else
            Axis = gca;
        end
    end
    
    ax = Axis;
    

end

