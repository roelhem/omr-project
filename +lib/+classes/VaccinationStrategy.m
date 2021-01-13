classdef VaccinationStrategy < lib.classes.AgeGroupPopulation
    %VACCINATIONSTRATEGY Holds the information of a vaccination stategy.
    %   Contains some methods to check if the VaccinationStrategy is
    %   allowed by the restrictions and methods to calculate the values
    %   that can be inserted into the SIRV_model.
    
    %% Initialisation
    methods
        function obj = VaccinationStrategy(scale, Rho, T)
            %VACCINATIONSTRATEGY Creates a new instance of this vaccination
            %   strategy.
            assert(isa(scale, 'lib.classes.AgeGroupPopulation'));
            obj = obj@lib.classes.AgeGroupPopulation(...
                scale.Boundaries, ...
                scale.Population ...
            );
            
            % Getting the default data.
            if nargin < 2
                Rho = zeros(obj.m, 0);
            end
            if nargin < 3
                T = 1;
            end
            
            if ~isduration(T)
                T = days(T);
            end
            
            % Initialize the timetable
            if isa(Rho, 'timetable')
                assert(width(Rho) == obj.m);
                obj.TTable = Rho;
            else
                assert(height(Rho) == obj.m);
                if width(Rho) == size(T)
                    obj.TTable = timetable( ...
                        'Size', [width(Rho) obj.m], ...
                        'VariableTypes', repmat("double", obj.m, 1), ...
                        'RowTimes', T ...
                    );
                elseif size(T) == 1
                    obj.TTable = timetable( ...
                        'Size', [width(Rho) obj.m], ...
                        'VariableTypes', repmat("double", obj.m, 1), ...
                        'TimeStep', T ...
                    );
                else
                    error('Invalid Time vector.');
                end
                
                obj.TTable{:,:} = Rho';
            end
            
            % Set the timetable variable names.
            obj.TTable.Properties.VariableNames = obj.GroupName;
            
            % Set the default InitialV.
            obj.InitialV = zeros(obj.m, 1);
        end
    end
    
    %% Stored properties
    properties
        TTable timetable % The underlying timetable holding the Vaccination
                         % Strategy data.
        InitialV double  % [m x 1] The initial amount of vaccinated people
                         %         per group j.
    end
    
    %% Dependent properties
    properties(Dependent)
        n      double    % [1 x 1] Number of timestamps for which a
                         %         value of the vaccination strategy is
                         %         known.
        Rho    double    % [m x n] Number of people that move from U to V
                         %         in group m at moment n due to this
                         %         vaccination strategy per day.
        RhoTot double    % [1 x n] The total of Rho over all groups.
        Nu     double    % [m x n] The nu that can be used in our model.
        P      double    % [m x n] The total amount of people that will
                         %         be vaccinated before the vaccination
                         %         moment.
        Ptot   double    % [m x 1] The total amount of people that will
                         %         be vaccinated at the and of this
                         %         vaccination strategy.
        T      duration  % [1 x n] The time after the initial state on
                         %         which the vaccination n will take
                         %         placed.
        DeltaT   double  % [1 x n] The amount of days from the previous
                         %         timestep, starting from 0.
        N        double  % [m x 1] Alias for the population of the age 
                         %         groups.
        InitialU double  % [m x 1] The initial amount of unvaccinated
                         %         people per group j.
        V        double  % [m x n] The amount of vaccinated people in
                         %         group j before datapoint n.
        U        double  % [m x n] The amount of unvaccinated people in
                         %         group j before datapoint n.
    end
    
    methods
        function out = get.n(obj)
            out = height(obj.TTable);
        end
        
        function out = get.Rho(obj)
            out = obj.TTable{:,:}';
        end
        
        function obj = set.Rho(obj, value)
            obj.TTable{:,:} = value';
        end
        
        function out = get.RhoTot(obj)
            out = sum(obj.Rho, 1);
        end
        
        function out = get.Nu(obj)
            out = obj.Rho ./ obj.U;
        end
        
        function out = get.P(obj)
            out = zeros(obj.m, obj.n);
            for i = 1:obj.n-1
                out(:,i+1) = out(:,i) + obj.Rho(:,i) .* obj.DeltaT(:,i);
            end
        end
        
        function obj = set.P(obj, value)
            assert(height(value) == obj.m);
            assert(width(value) == obj.n);
            for i = 1:obj.n-1
                obj.Rho(:,i) = (value(:,i+1) - value(:,i)) ./ obj.DeltaT(:,i);
            end
        end
        
        function out = get.Ptot(obj)
            out = sum(table2array(obj.DailyTable))';
        end
        
        function out = get.T(obj)
            out = obj.TTable.Properties.RowTimes';
        end
        
        function obj = set.T(obj, value)
            obj.TTable.Properties.RowTimes = value';
        end
        
        function out = get.DeltaT(obj)
            out = (obj.T(2:end) - obj.T(1:end-1)) / days(1);
        end
        
        function out = get.N(obj)
            out = obj.Population;
        end
        
        function obj = set.N(obj, value)
            assert(height(value) == obj.m);
            assert(width(value) == 1);
            obj.Population = value;
        end
        
        function out = get.InitialU(obj)
            out = obj.N - obj.InitialV;
        end
        
        function obj = set.InitialU(obj, value)
            assert(height(value) == obj.m);
            assert(width(value) == 1);
            obj.InitialV = obj.N - value;
        end
        
        function out = get.V(obj)
            out = repmat(obj.InitialV, 1, obj.n) + obj.P;
        end
        
        function obj = set.V(obj, value)
            obj.InitialV = value(:,1);
            obj.P = repmat(obj.InitialV, 1, obj.n) - value;
        end
        
        function out = get.U(obj)
            out = repmat(obj.N, 1, obj.n) - obj.V;
        end
        
        function obj = set.U(obj, value)
            obj.V = repmat(obj.N, 1, obj.n) - value;
        end
    end
    
    %% Plotting the vaccination strategy.
    methods
        function out = plotArea(obj)
            out = area(obj.T, obj.Rho' ...
            );
            colororder(jet(obj.m));
            ylabel('Amount of people per day');
            legend(obj.GroupCategory);
        end
        
        function out = plotLine(obj)
            out = plot(obj.T, obj.Rho');
        end
    end
    
    %% Mutating the vaccination strategy.
    methods
        
        function out = mutateDistMat(obj, grp, delta, add_dist, cmp_dist)
            out = lib.mutations.mut_dist_base(obj.Rho, grp, delta, add_dist, cmp_dist);
            
            % Check the output.
            valid_grp_total = all(abs(sum(out, 2) - sum(rho, 2)) < 1e-3);
            valid_day_total = all(abs(sum(out, 1) - sum(rho, 1)) < 1e-3);
            all_positive = all(out > -1e-5, 'all');
            
            if ~(valid_grp_total && valid_day_total && all_positive)
                disp('INVALID TRANSFORMATION!');
                
                disp('grp');
                disp(grp);
                
                disp('add_dist');
                disp(add_dist);
                
                disp('cmp_dist');
                disp(cmp_dist);
                
                disp('grp_diff:');
                disp(sum(out, 2) - sum(rho, 2));
                
                disp('day_diff:');
                disp(sum(out, 1) - sum(rho, 1));
                
                disp('delta:');
                disp([input_delta delta]);
                
                disp('compensate_add:');
                disp(compensate_add);
                
                disp('compensate_rem:');
                disp(compensate_rem);
                
                disp('out:');
                disp(out);
                
                disp('rho:');
                disp(rho);
                
%                  assignin('base', 'fcnState', struct( ...
%                      'vs', obj, ...
%                      'rho', rho, ...
%                      'out', out, ...
%                      'grp', grp, ...
%                      'add_dist', add_dist, ...
%                      'cmp_dist', cmp_dist, ...
%                      'input_add_dist', input_add_dist, ...
%                      'input_cmp_dist', input_cmp_dist, ...
%                      'grp_diff', sum(out, 2) - sum(rho, 2), ...
%                      'day_diff', sum(out, 2) - sum(rho, 2), ...
%                      'input_delta', input_delta, ...
%                      'adj_delta', delta, ...
%                      'compensate_add', compensate_add, ...
%                      'compensate_rem', compensate_rem, ...
%                      'valid_grp_total', valid_grp_total, ...
%                      'valid_day_total', valid_day_total, ...
%                      'all_positive', all_positive ...
%                  ));
                
                
                error('INVALID TRANSFORMATION!');
            end
        end
        
        function out = mutateSubstractMat(obj, grp, day, subs)
            out = obj.Rho;
            
            % Get the indices from which to compensate the mutation.
            comp_days = day+1:obj.n;
            comp_grps = 1:obj.m ~= grp;
            
            % Get total of the values after the mutated day.
            comp_days_total = sum(obj.Rho(:, comp_days), 2);
            comp_grps = and(comp_grps, comp_days_total' > 0);
            comp_total = sum(comp_days_total(comp_grps));
            
            % Compute the amount to remove from the provided group.
            subs = min([subs, obj.Rho(grp, day), comp_total]);
            
            if subs <= 1e-5
                return
            end
            
            % Compenstate for day total.
            delta_day = zeros(obj.m, 1);
            delta_day(grp) = -subs;
            delta_day(comp_grps, 1) = (comp_days_total(comp_grps) / comp_total) * subs;
            out(:, day) = out(:, day) + delta_day;
            
            % Compensate for group total.
            delta_comp_days = zeros(obj.m, length(comp_days));
            delta_comp_days(comp_grps, :) = - diag(delta_day(comp_grps)./comp_days_total(comp_grps)) * out(comp_grps, comp_days);
            delta_comp_days(grp, :) = - sum(delta_comp_days, 1);
            out(:, comp_days) = out(:, comp_days) + delta_comp_days;
        end
        
        function out = mutateAddMat(obj, grp, day, add)
            out = obj.Rho;
            
            % Get the indices from which to compensate the mutation.
            comp_days = day+1:obj.n;
            comp_grps = 1:obj.m ~= grp;
            
            % Get total of the values after the mutated day.
            comp_grps = and(comp_grps, obj.Rho(:, day)' > 0);
            comp_grps_day_total = sum(obj.Rho(comp_grps, day));
            mut_days_total = sum(obj.Rho(grp,comp_days));
            
            if comp_grps_day_total <= 1e-5
                return
            end
            
            % Compute the amount to remove from the provided group.
            add = min([add, comp_grps_day_total, mut_days_total]);
            
            if add <= 1e-5
                return
            end
            
            % Compenstate for day total.
            delta_day = zeros(obj.m, 1);
            delta_day(grp) = add;
            delta_day(comp_grps) = - (obj.Rho(comp_grps, day)) * add / comp_grps_day_total;
            out(:, day) = out(:, day) + delta_day;
            
            % Compensate for group total.
            delta_comp_days = zeros(obj.m, length(comp_days));
            delta_comp_days(grp, :) = obj.Rho(grp,comp_days) * add / mut_days_total;
            delta_comp_days(comp_grps, :) =  (delta_day(comp_grps) / add) * delta_comp_days(grp, :);
            out(:, comp_days) = out(:, comp_days) - delta_comp_days(:,:);
        end
        
        function out = mutateMat(obj, grp, day, delta)
            if delta > 0
                out = obj.mutateAddMat(grp, day, delta);
            elseif delta < 0
                out = obj.mutateSubstractMat(grp, day, -delta);
            else
                out = obj.Rho;
            end
        end
        
        function out = mutate(obj, grp, day, delta)
            out = obj;
            out.Rho = obj.mutateMat(grp, day, delta);
        end
        
        function out = mutateDist(obj, grp, delta, add_dist, cmp_dist)
            out = obj;
            out.Rho = obj.mutateDistMat(grp, delta, add_dist, cmp_dist);
        end
    end
    
    %% Reformatting the vaccination strategy.
    properties(Dependent)
        DailyTable duration % Timetable with a row for each day.
    end
    
    methods
        function out = get.DailyTable(obj)
            out = retime(obj.TTable, 'daily','mean');
            out = fillmissing(out, 'constant', 0);
        end
    end
    
    methods
        function result = V_at(obj, duration)
            % V_AT Returns the value of V at the provided duration.
            duration = days(duration);
            closestIndex = find(obj.T <= duration, 1, 'last');
            closestTime = obj.T(closestIndex);
            
            % Just give the InitialV if no times were available before
            % the provided time, just give the initial V.
            if(isempty(closestTime))
                result = obj.InitialV;
                return
            end
            
            % Compute the difference between the closest and the wanted
            % timestep.
            diffClosest = (duration) - closestTime;
            
            % If the difference is bigger than one day, it means that
            % we won't vaccinate on the days inbetween. Therefore, we can
            % just collapse the difference to one day.
            if(diffClosest >= days(1))
                diffClosest = days(1);
            end
            
            % Calculate the result.
            closestRho = obj.Rho(:,closestIndex);
            closestV   = obj.V(:,closestIndex);
            result = closestV + closestRho * diffClosest / days(1);
        end
        
        function result = retime(obj, NewT)
            % Get the new length of datapoints.
            if isduration(NewT)
                NewT = NewT/days(1);
            end
            New_n = length(NewT);
            
            % Compute the new Rho.
            NewRho = zeros(obj.m, New_n);
            for i = 1:New_n
                d = days(floor(NewT(i)));
                v = obj.DailyTable{d, :}';
                if ~isempty(v)
                    NewRho(:,i) = v;
                end
            end
            
            % Compute the new initial V.
            NewInitialV = obj.V_at(NewT(1));
            
            % Construct the new result.
            result = obj;
            result.InitialV = NewInitialV;
            result.TTable = array2timetable(NewRho', ...
                'RowTimes', days(NewT), ...
                'VariableNames', obj.GroupName ...
            );
        end
    end
    
    %% Convert strategy.
    methods
        function out = asInput(obj, n)
            if nargin < 2 || isempty(n)
                out = reshape(obj.Rho,obj.n * obj.m,1);
            else
                out = zeros(n*obj.m, 1);
                l = min(obj.n, n);
                out(1:l*obj.m) = reshape(obj.Rho(:, 1:l),l * obj.m,1);
            end
        end
    end
    
    methods(Static)
        function obj = fromInput(scale, x)
            m = scale.m;
            rho = reshape(x, m, length(x)/m);
            obj = lib.classes.VaccinationStrategy(scale, rho);
        end
    end
    
end

