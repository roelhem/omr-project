classdef AgeGroup
    %AGEGROUP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Boundaries
    end
    
    methods
        function obj = AgeGroup(input)
            %AGEGROUP Construct an instance of this class
            %   It can handle an array with boundaries or a string.
            if isnumeric(input)
                if width(input) == 1
                    obj.Boundaries = input;
                    obj.Boundaries(:,2) = Inf;
                else
                    obj.Boundaries = input;
                end
            else
                obj.Boundaries = lib.utils.strToBoundaries(input);
            end
        end
        
        function out = categories(obj)
            %METHOD Returns the categories of these age groups.
            out = lib.utils.boundariesToCat(obj.Boundaries);
        end
    end
end

