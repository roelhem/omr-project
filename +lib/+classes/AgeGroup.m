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
            if isa(input, 'lib.classes.AgeGroup')
                obj.Boundaries = input.Boundaries;
            elseif isnumeric(input)
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
            %CATEGORIES Returns the categories of these age groups.
            out = lib.utils.boundariesToCat(obj.Boundaries);
        end
        
        function out = table(obj)
            %TABLE Returns the data of this AgeGroup class as a table.
            Categories = obj.categories;
            out = table(Categories, obj.Boundaries,...
                'VariableNames', {'Category', 'Boundaries'},...
                'RowNames', strings(Categories)...
            );
        end
        
        function out = size(obj)
            %SIZE Returns the amount of age groups stored.
            out = height(obj.Boundaries);
        end
        
        function out = groupOverlap(obj, V)
            %GROUPOVERLAP Returns the proportion of overlap of a certain
            %group with a group managed in this class.
            
            % Ensure that the dimensions of the provided groups are right.
            assert(width(V) == 2, "Given valueGroupBoundaries doesn't have width 2.")
            
            % Init the result.
            m = height(V);
            out = zeros(obj.size, m);
            
            % Set infinite values to 100.
            B = obj.Boundaries;
            B(isinf(B)) = 99;
            V(isinf(V)) = 99;
            
            % Make upper bound open.
            B(:,2) = B(:,2);
            V(:,2) = V(:,2);
            
            % Loop through all own groups B.
            for i = 1:obj.size                
                % Loop through all given groups V.
                for j = 1:m
                    % Get the size of the overlap.
                    maxLow  = max(B(i,1), V(j,1)); % The maximum upper bound.
                    minUp   = min(B(i,2), V(j,2)); % The minimum lower bound.
                    overlap = minUp - maxLow + 1; % The size of the overlap.
                    % Check if B and V overlap.
                    if overlap > 0
                        % Set the result value to the proportion of V that
                        % overlaps with B.
                        Vsize = V(j,2) - V(j,1) + 1;
                        out(i,j) = overlap / Vsize;
                    end
                end
            end
            
        end
        
        function out = sumPerGroup(obj, values, valueGroupBoundaries)
            %SUMPERGROUP Takes a vector of values and the corresponding
            %   boundaries of the groups. It will sum the values per 
            %   managed age group.
            
            % Ensure that the input has the right dimensions.
            B = lib.classes.AgeGroup(valueGroupBoundaries);
            assert(height(values) == height(B.Boundaries), "Height of values and valueGroupBoundaries don't match.")
            
            % Initialize the result.
            out = obj.groupOverlap(B.Boundaries) * values;
        end
        
        function out = avgResize(obj, values, valueGroupBoundaries)
            %AVGRESIZE Takes a vector of values and resizes the values
            %   by taking the weighted average of the overlap.
            
            % Ensure that the input has the right dimensions.
            B = lib.classes.AgeGroup(valueGroupBoundaries);
            assert(height(values) == height(B.Boundaries), "Height of values and valueGroupBoundaries don't match.")
            
            % Get the overlap matrix.
            Overlap = obj.groupOverlap(B.Boundaries);
            % Initialize the result.
            out = (Overlap * values) ./ sum(Overlap, 2);
        end
        
        function out = resize(obj, V, B)
            out = obj.sumPerGroup(V, B);
        end
    end
end

