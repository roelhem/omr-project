classdef ModelState < lib.classes.AgeGroupPopulation
    %MODELSTATE Manages all the variables of a model at some timestap.
    
    %% Initialisation methods
    methods
        function obj = ModelState(group, N)
            %MODELSTATE Initialize a new ModelState object.
            
            % Initialize the superclass.
            obj = obj@lib.classes.AgeGroupPopulation(group, N);
            
            % Set default values.
            obj.T         = days(0);
            obj.S         = obj.N;
            obj.I         = zeros(obj.m, 1);
            obj.R         = zeros(obj.m, 1);
            obj.Alpha     = ones(obj.m, 1) * 0.1;
            obj.Beta      = eye(obj.m, obj.m);
            obj.Nu        = zeros(obj.m, 1);
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
            obj = obj.loadInitialValues(D)...
                    .loadContactMatrix(CFile)...
                    .loadReprNumData(D);
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
        T         % [1 x 1] The amount of days from the initial state
                  %         of the model.
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
        N          % [m x 1] The total amount of people per age group.
        V          % [m x 1] The number of vaccinated people per age group.
        U          % [m x 1] The number of unvaccinated people per age.
                   % group.
        Tau        % [m x 1] The average length of time for which an infected 
                   % individual remains infectious to others.
        Rho        % [m x 1] The amount of people that we are planning to
                   %         vaccinate in group j per day.
        C          % [m x m] The contact matrix, which is the chance that
        P_CtoI     % [m x m] The probability that a contact between an infected
                   % individual and a susceptible individual will lead to a
                   % new infection.
        ReprNum    % [1 x 1] The average number of people that an infectious person
                   % can be expected to infect before they move to the
                   % 'removed' compartment, assuming that all their contacts
                   % are with people in the susceptible compartment.
        ReprEff    % [1 x 1] The effective reproduction number.
        AsympRatio % [m x 1] The ratio of people that are infectious, but
                   % don't have symptoms.
        I_Symp     % [m x 1] The number of people that are infectious and
                   % have symptoms.
        I_Asymp    % [m x 1] The number of people that are infectious and
                   % don't have symptoms.
                   
        StoI       % [m x 1] The amount of people that move from S to I in
                   % one day per age group.
        ItoR       % [m x 1] The amount of people that move from I to R in
                   % one day per age group.
        StoV       % [m x 1] The amount of people that move from S to V in
                   % one day per age group.
        ItoV       % [m x 1] The amount of people that move from I to V in
                   % one day per age group.
        RtoV       % [m x 1] The amount of people that move from R to V in
                   % one day per age group.
        UtoV       % [m x 1] The amount of unvaccinated people that will
                   % be vaccinated in one day per age group.
                   
        dS         % [m x 1] The change of S in people per day.
        dI         % [m x 1] The change of I in people per day.
        dR         % [m x 1] The change of R in people per day.
        dU         % [m x 1] The change of U in people per day.
        dV         % [m x 1] The change of V in people per day.
    end
    
    methods
        
        function out = get.N(obj)
            % GET.N Alias for the population per age group.
            out = obj.Population;
        end
        
        function obj = set.N(obj, value)
            % SET.N Setting the population per age group.
            %       This is just an alias for the `Population` on the
            %       underlying age group distribution.
            assert(height(value) == obj.m);
            obj.Population = value;
        end
        
        function out = get.U(obj)
            % GET.U Gets the total amount of unvaccinated people.
            %       Can be retrieved by getting the sum of `S`, `I` and
            %       `R`.
            out = obj.S + obj.I + obj.R;
        end
        
        function obj = set.U(obj, value)
            % SET.U Sets the total amount of unvaccinated people.
            %       It will adjust `S`, `I` and `R` so that the ratios
            %       Between those values will stay the same.
            assert(height(value) == obj.m);
            PrevU = obj.U;
            DeltaU = value - PrevU;
            obj.S = obj.S + DeltaU * obj.S ./ PrevU;
            obj.I = obj.I + DeltaU * obj.I ./ PrevU;
            obj.R = obj.R + DeltaU * obj.R ./ PrevU;
        end
        
        function out = get.V(obj)
            % GET.V Gets the total amount of vaccinated people per age
            %       group.
            %       It will calculate this by substracting `U` from `N`.
            out = obj.N - obj.U;
        end
        
        function obj = set.V(obj, value)
            % SET.V Sets the total amount of vaccinated people per age
            %       group. It will calculate this by setting `U`.
            %       @see `set.U`.
            assert(height(value) == obj.m);
            obj.U = obj.N - value;
        end
        
        function out = get.Tau(obj)
            out = 1 ./ obj.Alpha;
        end
        
        function obj = set.Tau(obj, value)
            if(size(value) == 1)
                value = ones(obj.m, 1) * value;
            elseif(height(value) ~= obj.value || width(value) ~= 1)
                error('Tau has wrong size.');
            end
            obj.Alpha = 1 ./ value;
        end
        
        function out = get.Rho(obj)
            out = obj.Nu .* obj.U;
        end
        
        function obj = set.Rho(obj, value)
            obj.Nu = value ./  obj.U;
        end
        
        function out = get.C(obj)
            out = obj.Beta * diag(1 ./ obj.SusVect);
        end
        
        function obj = set.C(obj, value)
            assert(width(value) == obj.m);
            assert(height(value) == obj.m);
            obj.Beta =  diag(obj.SusVect) * value;
        end
        
        function out = get.P_CtoI(obj)
            out = obj.Beta ./ obj.Alpha; % TODO: Check this.
        end
        
        function obj = set.P_CtoI(obj, value)
            obj.Beta = value .* obj.Alpha; % TODO: Check this.
        end
        
        function out = get.ReprNum(obj)
            % GET.REPRNUM Gets the reproduction number, when everyone is
            %             is susceptible.
            %             The value can is based on the identity.
            %             ```LaTeX
            %             \R_0 = \frac{\beta}{\alpha} \cdot N.
            %             ```
            out = max(eig(obj.Beta)) / mean(obj.Alpha);
        end
        
        function obj = set.ReprNum(obj, value)
            % SET.REPRNUM Ajusts `Beta` such that the reproduction number,
            %             will be the given value.
            %             The value can is based on the identity.
            %             ```LaTeX
            %             \R_0 = \frac{\beta}{\alpha} \cdot N.
            %             ```
            assert(width(value) == 1);
            assert(height(value) == 1);
            obj.Beta = obj.Beta * (value / obj.ReprNum);
        end
        
        function out = get.ReprEff(obj)
            % GET.REPREFF Gets the effective reproduction number `R_eff`.
            %             The value can is based on the identity.
            %             ```LaTeX
            %             \R_eff = \frac{\beta}{\alpha} \cdot S.
            %             ```
            out = max(eig(obj.Beta)) / mean(obj.Alpha) * (sum(obj.S) / sum(obj.N));
        end
        
        function obj = set.ReprEff(obj, value)
            % SET.REPREFF Ajusts the `Beta` so that the effective
            %             reproduction number becomes `R_eff`.
            assert(width(value) == 1);
            assert(height(value) == 1);
            obj.Beta = obj.Beta * (value / obj.ReprEff);
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
        
        function out = get.StoI(obj)
            out = (obj.Beta * (obj.I ./ obj.N)) .* obj.S;
        end
        
        function out = get.ItoR(obj)
            out = obj.Alpha .* obj.I;
        end
        
        function obj = set.ItoR(obj, value)
            assert(width(value) == 1);
            assert(height(value) == obj.m);
            obj.Alpha = value ./ obj.I;
        end
        
        function out = get.StoV(obj)
            out = obj.Nu .* obj.S;
        end
        
        function obj = set.StoV(obj, value)
            assert(width(value) == 1);
            assert(height(value) == obj.m);
            obj.Nu = value ./ obj.S;
        end
        
        function out = get.ItoV(obj)
            out = obj.Nu .* obj.I;
        end
        
        function obj = set.ItoV(obj, value)
            assert(width(value) == 1);
            assert(height(value) == obj.m);
            obj.Nu = value ./ obj.I;
        end
        
        function out = get.RtoV(obj)
            out = obj.Nu .* obj.R;
        end
        
        function obj = set.RtoV(obj, value)
            assert(width(value) == 1);
            assert(height(value) == obj.m);
            obj.Nu = value ./ obj.R;
        end
        
        function out = get.UtoV(obj)
            out = obj.Nu .* obj.U;
        end
        
        function obj = set.UtoV(obj, value)
            assert(width(value) == 1);
            assert(height(value) == obj.m);
            obj.Nu = value ./ obj.U;
        end
        
        function out = get.dS(obj)
            out = - obj.StoI - obj.StoV;
        end
        
        function out = get.dI(obj)
            out = obj.StoI - obj.ItoR - obj.ItoV;
        end
        
        function out = get.dR(obj)
            out = obj.ItoR - obj.RtoV;
        end
        
        function out = get.dU(obj)
            out = -obj.UtoV;
        end
        
        function out = get.dV(obj)
            out = obj.UtoV;
        end
    end
    
    %% Function handles.
    properties(Dependent)
        vacc_vectorized_run function_handle
        vacc_run function_handle
    end
    
    methods
        function M = getModelResult(obj, S_res, I_res, U_res, Nu_res)
            % Get the sizes.
            n = size(S_res,3);
            R_res = U_res - (S_res + I_res);

            % Get the states.
            State(1:n) = obj;
            for j = 1:n
                State(j).T(:) = obj.T + days(j - 1);
                State(j).S(:,:) = S_res(1,:,j);
                State(j).I(:,:) = I_res(1,:,j);
                State(j).R(:,:) = R_res(1,:,j);
                State(j).Nu(:,:) = Nu_res(1,:,j);
            end
            
            % Get the results.
            M = lib.classes.ModelResult(State, 1, "EulerForward");
        end
        
        function out = get.vacc_vectorized_run(obj)
            vectorized_fcn = lib.models.SIRV_EulerForward_vectorized(obj);
            
            function result = f(x)
                [S_res,I_res,U_res,StoI_res,Nu_res] = vectorized_fcn(x);
                
                p = size(S_res, 1);
                
                result = struct( ...
                    'p', p, ...
                    'S', S_res, ...
                    'I', I_res, ...
                    'U', U_res, ...
                    'StoI', StoI_res, ...
                    'Nu', Nu_res, ...
                    'getModelResult', @(i)obj.getModelResult(S_res(i,:,:), I_res(i,:,:), U_res(i,:,:),Nu_res(i,:,:)) ...
                );
            end
            
            out = @f;
            
        end
        
        function out = get.vacc_run(obj)
            aa = lib.models.SIRV_EulerForward_vectorized(obj);
            out = @(x)aa(x');
        end
    end
    
    
    %% Loading from datasources.
    methods
        function obj = loadVaccinationStrategy(obj, Strategy, TimeStep)
            % Initialize the strategy if nessicary.
            if ~isa(Strategy, 'lib.classes.VaccinationStrategy')
                if nargin < 3
                    TimeStep = days(1);
                end
                Strategy = lib.classes.VaccinationStrategy(obj, Stategy, TimeStep);
            end
            
            % Retime and set strategy to the time of this ModelState.
            v = Strategy.retime(obj.T);
            v.InitialV = obj.V;
            
            % Set the Rho of this ModelState.
            obj.Rho = v.Rho(:,1);
        end
        
        function obj = loadContactMatrix(obj, filename)
            obj.C = lib.conversions.contactMatrix(...
                filename,...
                obj.Boundaries...
            );
        end
        
        function obj = loadReprNumData(obj, at_date)
            global rivm_reproduction
            
            Filtered = rmmissing(rivm_reproduction);
            Filtered = Filtered(Filtered.Date <= datetime(at_date), :);
            obj.ReprEff = Filtered{end, 'Rt_avg'};
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
    
    %% Variable conversions.
    methods
        function out = toVaccinationStrategy(obj)
            out = lib.classes.VaccinationStrategy( ...
                obj, ...
                obj.Rho, ...
                obj.T ...
            );
            out.InitialV = obj.V;
        end
        
        function out = toModelResult(obj, DeltaT, Method)
            % Get the default variables.
            if nargin <= 2
                DeltaT = 1;
            end
            if nargin <= 3
                Method = "EulerForward";
            end
            
            % Return the ModelResult.
            out = lib.classes.ModelResult( ...
                obj, ...
                DeltaT, ...
                Method ...
            );
        end
    end
    
    %% Running the model.
    methods
        
        function result = run(obj, DeltaT, n, Method, VaccinationStrategy)
            % RUN Runs the model and return a `ModelResult instance`.
            
            % Set default argument values.
            if nargin <= 3
                Method = "EulerForward";
            end
            if nargin <= 4
                VaccinationStrategy = obj.toVaccinationStrategy();
            end
            
            % Check input arguments.
            assert(isa(VaccinationStrategy, 'lib.classes.VaccinationStrategy'), 'Invalid Vaccination strategy.');
            assert(DeltaT > 0, 'DeltaT has to be bigger than 0');
            assert(n > 0, 'n has to be bigger than 0');
            
            % Get Nu
            dt_days = days(DeltaT);
            T_end = obj.T + dt_days * (n - 1);
            nu = VaccinationStrategy.retime(obj.T:dt_days:T_end).Nu;
            
            % Run the model with the chosen method.
            switch string(Method)
                case "EulerForward"
                    [SRes, IRes, RRes] = obj.run_EulerForward(DeltaT, n, nu);
                case "EulerBackward"
                    [SRes, IRes, RRes] = obj.run_EulerBackward(DeltaT, n, nu);
                otherwise
                    error(strcat("Method '", string(Method), "' not recognised/implemented."));
            end
            
            % Create the state vector.
            State(1:n) = obj;
            for i = 1:n
                State(i).T  = obj.T + days(DeltaT * (i - 1));
                State(i).S  = SRes(:,i);
                State(i).I  = IRes(:,i);
                State(i).R  = RRes(:,i);
                State(i).Nu = nu(:,i);
            end
            
           
            % Create the model result.
            result = lib.classes.ModelResult(State, DeltaT, Method);
        end
        
        function result = next(obj, DeltaT, Method)
            % NEXT Get the next state with the provided DeltaT and Method.
            if nargin <= 1
                DeltaT = 1;
            end
            if nargin <= 2
                Method = "EulerForward";
            end
            
            assert(DeltaT > 0, 'DeltaT has to be bigger than 0');
            
            % Run the model with the chosen method.
            switch string(Method)
                case "EulerForward"
                    [SRes, IRes, RRes] = obj.run_EulerForward(DeltaT, 2);
                case "EulerBackward"
                    [SRes, IRes, RRes] = obj.run_EulerBackward(DeltaT, 2);
                otherwise
                    error(strcat("Method '", string(Method), "' not recognised/implemented."));
            end
            
            result = obj;
            result.T = obj.T + days(DeltaT);
            result.S = SRes(:,2);
            result.I = IRes(:,2);
            result.R = RRes(:,2);
        end
        
        function [S,I,R,V] = run_EulerForward(obj, DeltaT, n, nu)
            % RUN_EULERFORWARD Runs the model using the EulerForward
            %                  numeric method.
            if nargin < 4 || isempty(nu)
                nu = [obj.Nu, zeros(obj.m, n - 1)];
            end
            
            [SRes,IRes,RRes,VRes] = lib.models.SIRV_EulerForward(...
                obj.S,...
                obj.I,...
                obj.R,...
                obj.V,...
                obj.Beta,...
                obj.Alpha,...
                nu,...
                DeltaT,...
                n...
            );
        
            S = SRes;
            I = IRes;
            R = RRes;
            V = VRes;
        end
        
        function [S,I,R,V] = run_EulerBackward(obj, DeltaT, n, nu)
            % RUN_EULERBACKWARD Runs the model using the EulerBackward
            %                   numeric method.
            if nargin < 4 || isempty(nu)
                nu = [obj.Nu, zeros(obj.m, n - 1)];
            end
            
            [SRes,IRes,RRes,VRes] = lib.models.SIRV_EulerBackward(...
                obj.S,...
                obj.I,...
                obj.R,...
                obj.V,...
                obj.Beta,...
                obj.Alpha,...
                nu,...
                DeltaT,...
                n...
            );
        
            S = SRes;
            I = IRes;
            R = RRes;
            V = VRes;
        end
    end
    
    %% Helpers to build an UI.
    methods
        function T = buildGoupTable(obj)
            T = table( ...
                round(obj.N), ...
                round(obj.S), ...
                round(obj.I), ...
                round(obj.R), ...
                round(obj.V), ...
                obj.Alpha, ...
                obj.Tau, ...
                obj.SusVect, ...
                obj.SympRatio, ...
                round(obj.I_Symp), ...
                'VariableNames', [
                    "N_0"
                    "S_0"
                    "I_0"
                    "R_0"
                    "V_0"
                    "Alpha"
                    "Tau"
                    "SusVect"
                    "SympRatio"
                    "I_Symp_0"
                ], ...
                'RowNames', string(obj.categories()) ...
            );
        end
    end
    
    %% Plots
    methods
        function h = plotHeatMap(obj, varargin)
            par = gcf;
            
            Data = "Beta"; % "C" "P_CtoI"
            Title = [];
            XLabel = "Transmitting Group";
            YLabel = "Receiving Group";
            
            for ii = 1:2:length(varargin)
                switch string(varargin{ii})
                    case "Data"
                        Data = varargin{ii + 1};
                    case "XLabel"
                        XLabel = varargin{ii + 1};
                    case "YLabel"
                        YLabel = varargin{ii + 1};
                end
            end
            
            values = obj.GroupName;
            
            % Getting the heatmap.
            switch string(Data)
                case "Beta"
                    h = heatmap(par, values, values, obj.Beta);
                case "C"
                    h = heatmap(par, values, values, obj.C);
                case "P_CtoI"
                    h = heatmap(par, values, values, obj.P_CtoI);
                otherwise
                    error(strcat("Doesn't recognise data '", Data, "'."));
            end
            
            if ~isempty(Title)
                h.Title = Title;
            else
                switch string(Data)
                    case "Beta"
                        h.Title = "Values of \beta";
                    case "C"
                        h.Title = "Contact matrix (Adjusted for R_0)";
                    case "P_CtoI"
                        h.Title = "P_{C \to I}";
                end
            end
            
            h.XLabel = XLabel;
            h.YLabel = YLabel;
            
        end
    end
end

