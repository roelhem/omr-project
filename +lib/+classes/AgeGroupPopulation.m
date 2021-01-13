classdef AgeGroupPopulation < lib.classes.AgeGroup
    %AGEGROUPPOPULATION Manages age groups and there associate population.
    
    properties
        Population
    end
    
    methods(Static)
        function obj = fromCbs(a, b)
            global cbs_populationTotal;
            
            switch nargin
                case 2
                    groupCategory = a;
                    gender = b;
                case 1
                    groupCategory = a;
                    gender = 'Total';
                otherwise
                    groupCategory = 'One';
                    gender = 'Total';
            end
            
            T = cbs_populationTotal.(groupCategory);
            obj = lib.classes.AgeGroupPopulation(...
                T.Group.Boundaries,...
                T.(gender)...
            );
        end
    end
    
    methods
        function obj = AgeGroupPopulation(groups, population)
            %AGEGROUPPOPULATION Construct an instance of this class
            
            if nargin == 0
                groups = "0+";
                population = 1;
            end
            
            obj = obj@lib.classes.AgeGroup(groups);
            
            assert(width(population) == 1, "Population must have width 1");
            assert(height(population) == obj.m, "Population must have the same dimensions as the amount of groups.");
            assert(isnumeric(population), "Population must be an array of numeric values.")
            
            obj.Population = population;
        end
        
        function out = rescale(obj, groups)
            %RESCALE Scales these age-groups according to the provided
            %groups.
            other = lib.classes.AgeGroup(groups);
            out = lib.classes.AgeGroupPopulation(...
                groups,...
                other.sumPerGroup(obj.Population, obj.Boundaries)...
            );
        end
        
        function out = weightedGroupOverlap(obj, G)
            %WEIGHTEDGROUPOVERLAP Determines a matrix that transforms
            %   data with other group distributions to group distributions
            %   that are compatible with the groups that this
            %   AgeGroupPopulation instance manages.
            
            % Ensure that the input is a AgeGroupPopulation instance.
            assert(isa(G, 'lib.classes.AgeGroupPopulation'), "Groups must be an instance of lib.classes.AgeGroupPopulation");
            
            % Get the overlap matrix.
            O = obj.groupOverlap(G.Boundaries);
            
            % Get the matrix that weights the values with the population.
            % out = diag(1./obj.Population) * O * diag(G.Population);
            out = diag(1./obj.Population) * O * diag(G.Population);
        end
        
        function out = weightedResize(obj, V, G)
            %WEIGHTEDRESIZE Weights
            
            % Initialize the result.
            O = obj.weightedGroupOverlap(G);
            out = O * V;
        end
    end
    
    %% Dependent properties
    properties(Dependent)
        PopulationTotal double
    end
    
    methods
        function out = get.PopulationTotal(obj)
            out = sum(obj.Population);
        end
    end
    
    
end

