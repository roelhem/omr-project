classdef ModelState < lib.classes.AgeGroupPopulation
    %MODELSTATE Manages all the variables of a model at some timestap.
    
    %% Initialisation methods
    methods
        function obj = ModelState(group, N)
            %MODELSTATE Initialize a new ModelState object.
            
            % Initialize the superclass.
            obj = obj@lib.classes.AgeGroupPopulation(group, N);
            
            % Set default values.
            obj.S         = obj.N;
            obj.I         = zeros(obj.m, 1);
            obj.R         = zeros(obj.m, 1);
            obj.Alpha     = ones(obj.m, 1) * 0.1;
            obj.Beta      = eye(obj.m, obj.m);
            obj.Nu        = zeros(obj.m, 0);
            obj.SusVect   = obj.avgResize(obj.STD_SusVect, obj.STD_Boundaries);
            obj.SympRatio = obj.avgResize(obj.STD_SympRatio, obj.STD_Boundaries);
        end
    end
    
    methods(Static)
        function obj = standardData(at_date, groups, contactMatrix)
            % STANDARDDATA Creates a new ModelState with standard data.
            %
            %   Arguments (All optional.)
            %   `at_date`: The date from which you want to start the data.
            %              Defaults to one month ago.
            %   `groups`:  The groups in which you want to seperate the
            %              population.
            %              Defaults to `STD_Boundaries`.
            %   `contactMatrix`: The contact matrix that you want to load.
            %                    Defaults to 'Contact_matrix.csv'.
            
            global cbs_AgeGroupPopulation;
            
            % Parsing the arguments.
            D = datetime() - calmonths(1);
            G = lib.classes.ModelState.STD_Boundaries;
            CFile = 'Contact_matrix.csv';
            if nargin >= 1
                D = datetime(at_date);
            end
            if nargin >= 2
                G = groups;
            end
            if nargin >= 2
                CFile = contactMatrix;
            end
            
            % Scale the population.
            scale = cbs_AgeGroupPopulation.rescale(G);
            
            % Create the object.
            obj = lib.classes.ModelState(scale.Boundaries, scale.Population);
            
            % Fill the object with data from the api's.
            obj = obj.loadContactMatrix(CFile)...
                    .loadReprNumData(D)...
                    .loadInitialValues(D);
        end
    end
    
    %% Standard Data.
    properties(Constant)
        STD_Boundaries = lib.utils.createGroupBoundaries(10, 70);
        STD_SusVect    = [0.4  ;0.38 ;0.79 ;0.86 ;0.8  ;0.82 ;0.88 ;0.74 ];
        STD_SympRatio  = [0.701;0.701;0.626;0.596;0.573;0.599;0.616;0.687];
    end
    
    %% Stateful properties.
    properties
        S         % [m x 1] The number of susseptible people per age group.
        I         % [m x 1] The number of infectious people per age group.
        R         % [m x 1] The number of recovered people per age group.
        Alpha     % [m x 1] The probability that an individual from the infectious 
                  % group moves to the removed compartment in a given 
                  % time-step.
        Beta      % [m x m] The average number of people that an infectious person
                  % in group j can be expected to infect in group i in a
                  % assuming that all their contacts are with people in the 
                  % susceptible compartment.
        Nu        % [m x 1] The proportion unvaccinated people that we want to vaccinate
                  % in the current timestep.
        SusVect   % [m x 1] The susceptibility of a person in group i.
        SympRatio % [m x 1] The number of people that actually have symptoms when
                  % they are infectious.
    end
    
    %% Derived properties.
    properties (Dependent)
        m          % The amount of age groups.
        N          % [m x 1] The total amount of people per age group.
        V          % [m x 1] The number of vaccinated people per age group.
        U          % [m x 1] The number of unvaccinated people per age.
                   % group.
        Tau        % [m x 1] The average length of time for which an infected 
                   % individual remains infectious to others.
        C          % [m x m] The contact matrix, which is the chance that
        P_CtoI     % [m x m] The probability that a contact between an infected
                   % individual and a susceptible individual will lead to a
                   % new infection.
        ReprNum    % [1 x 1] The average number of people that an infectious person
                   % can be expected to infect before they move to the
                   % 'removed' compartment, assuming that all their contacts
                   % are with people in the susceptible compartment.
        AsympRatio % [m x 1] The ratio of people that are infectious, but
                   % don't have symptoms.
        I_Symp     % [m x 1] The number of people that are infectious and
                   % have symptoms.
        I_Asymp    % [m x 1] The number of people that are infectious and
                   % don't have symptoms.
    end
    
    methods
        function out = get.m(obj)
            out = obj.size;
        end
        
        function out = get.N(obj)
            out = obj.Population;
        end
        
        function obj = set.N(obj, value)
            assert(height(value) == obj.m);
            obj.population = value;
        end
        
        function out = get.U(obj)
            out = obj.S + obj.I + obj.R;
        end
        
        function obj = set.U(obj, value)
            assert(height(value) == obj.m);
            PrevU = obj.U;
            DeltaU = value - PrevU;
            obj.S = obj.S + DeltaU * obj.S ./ PrevU;
            obj.I = obj.I + DeltaU * obj.I ./ PrevU;
            obj.R = obj.R + DeltaU * obj.R ./ PrevU;
        end
        
        function out = get.V(obj)
            out = obj.N - obj.U;
        end
        
        function obj = set.V(obj, value)
            assert(height(value) == obj.m);
            obj.U = obj.N - value;
        end
        
        function out = get.Tau(obj)
            out = 1 ./ obj.Alpha;
        end
        
        function obj = set.Tau(obj, value)
            if(size(value) == 1)
                value = ones(obj.size, 1) * value;
            elseif(height(value) ~= obj.value || width(value) ~= 1)
                error('Tau has wrong size.');
            end
            obj.Alpha = 1 ./ value;
        end
        
        function out = get.C(obj)
            out = diag(1 ./ obj.SusVect) * obj.P_CtoI;
        end
        
        function obj = set.C(obj, value)
            assert(width(value) == obj.m);
            assert(height(value) == obj.m);
            obj.Beta = value * diag(obj.SusVect);
        end
        
        function out = get.P_CtoI(obj)
            out = obj.Beta; % TODO: Check this.
        end
        
        function obj = set.P_CtoI(obj, value)
            obj.Beta = value; % TODO: Check this.
        end
        
        function out = get.ReprNum(obj)
            out = max(eig(obj.Beta));
        end
        
        function obj = set.ReprNum(obj, value)
            assert(width(value) == 1);
            assert(height(value) == 1);
            obj.Beta = obj.Beta * (value / obj.ReprNum);
        end
        
        function out = get.AsympRatio(obj)
            out = ones(obj.m, 1) - obj.SympRatio;
        end
        
        function obj = set.AsympRatio(obj, value)
            assert(width(value) == 1);
            assert(height(value) == obj.m);
            assert(all(value <= 1) && all(value >= 0));
            obj.SympRatio = ones(obj.m, 1) - value;
        end
        
        function out = get.I_Symp(obj)
            out = obj.I .* obj.SympRatio;
        end
        
        function obj = set.I_Symp(obj, value)
            assert(width(value) == 1);
            assert(height(value) == obj.m);
            obj.I = value ./ obj.SympRatio;
        end
        
        function out = get.I_Asymp(obj)
            out = obj.I .* obj.AsympRatio;
        end
        
        function obj = set.I_Asymp(obj, value)
            assert(width(value) == 1);
            assert(height(value) == obj.m);
            obj.I = value ./ obj.AsympRatio;
        end
    end
    
    %% Loading from datasources.
    methods
        function obj = loadContactMatrix(obj, filename)
            obj.C = lib.conversions.contactMatrix(...
                filename,...
                obj.Boundaries...
            );
        end
        
        function obj = loadReprNumData(obj, at_date)
            global rivm_reproduction
            
            T = rmmissing(rivm_reproduction);
            T = T(T.Date <= datetime(at_date), :);
            obj.ReprNum = T{end, 'Rt_avg'};
        end
        
        function obj = loadInitialValues(obj, at_date)
            global cbs_AgeGroupPopulation;
            
            exp_R_date = datetime(at_date) - days(mean(obj.Tau));
            
            removed = lib.conversions.casesPerGroup(exp_R_date);
            infectious = lib.conversions.casesPerGroup(exp_R_date, at_date);
            
            % Get the age group scales of the casesPerGroup data for
            % rescaling.
            Rscale = cbs_AgeGroupPopulation.rescale(removed.Agegroup);
            if removed.Agegroup == infectious.Agegroup
                % Just reuse the Rscale if the scales match.
                Iscale = Rscale;
            else
                Iscale = cbs_AgeGroupPopulation.rescale(infectious.Agegroup);
            end
            
            % Resize the groups and set the values.
            obj.R = obj.weightedResize(removed.GroupCount, Rscale);
            obj.I_Asymp = obj.weightedResize(infectious.GroupCount, Iscale);
            obj.S = obj.N - (obj.R + obj.I);
        end
    end
    
    %% Running the model.
    methods
        function [S,I,R,V] = run_EulerForward(obj, DeltaT, n)
            [SRes,IRes,RRes,VRes] = lib.models.SIRV_EulerForward(...
                obj.S,...
                obj.I,...
                obj.R,...
                obj.V,...
                obj.Beta,...
                obj.Alpha,...
                zeros(obj.m, n),...
                DeltaT,...
                n...
            );
        
            S = SRes;
            I = IRes;
            R = RRes;
            V = VRes;
        end
    end
end

