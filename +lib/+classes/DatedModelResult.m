classdef DatedModelResult < lib.classes.ModelResult
    %DATEDMODELRESULT A ModelResult that knows from which date the result
    %   was running.
    
    %% Serialisation
    %MODELRESULT The results got form a model that run successfully.
    methods
        function out = toStruct(obj)
            out = obj.toStruct@lib.classes.ModelResult();
            out.StartDate = obj.StartDate;
        end
    end
    
    methods(Static)
        function out = fromStruct(s)
            State = arrayfun(@(x)lib.classes.ModelState.fromStruct(x), s.State);
            out = lib.classes.DatedModelResult(State, s.DeltaT, s.Method, s.StartDate);
        end
    end
    
    %% Initialisation
    properties
        StartDate datetime % [1 x 1] The date from which the model will 
                           %         start.
    end
    
    methods
        function obj = DatedModelResult(State, DeltaT, Method, StartDate)
            %DATEDMODELRESULT Creates a new Dated Model Result.
            obj = obj@lib.classes.ModelResult(State, DeltaT, Method);
            obj.StartDate = StartDate;
        end
    end
    
    %% Dependent properties
    properties(Dependent)
        Date datetime % [1 x n] The datetime that belongs to timestamp j.
    end
    
    methods
        function out = get.Date(obj)
            out = obj.StartDate + obj.T;
        end
    end
    
    %% Superclass overrides.
    methods
        
        function XData = getXData(obj, varargin)
            XData = 'Date';
            
            for ii = 1:2:length(varargin)
                switch string(varargin{ii})
                    case "XData"
                        XData = varargin{ii + 1};
                end
            end
            
            if isequal(XData, 'Date')
                XData = obj.Date;
            elseif isequal(XData, 'Time')
                XData = obj.T;
            end
            
        end
        
        function plotRepr(obj, varargin)
            XData = obj.getXData(varargin{:});
            
            ax = lib.plots.get_axis(varargin{:});
            ylabel(ax, 'Infections per infected individual');
            hold(ax, 'on');
            plot(ax, XData, obj.ReprNum, 'DisplayName', 'R_0');
            plot(ax, XData, obj.ReprEff, 'DisplayName', 'R_{eff}');
            legend(ax, 'show');
            title(ax, 'Reproduction rate');
            hold(ax, 'off');
        end
        
        function plotPerGroupLines(obj, varargin)
            Plots = {'S', 'I', 'R'};
            for ii = 1:2:length(varargin)
                switch string(varargin{ii})
                    case "Plots"
                        Plots = varargin{ii + 1};
                end
            end
            
            XData = obj.getXData(varargin{:});
            
            ax = lib.plots.get_axis(varargin{:});
            ylabel(ax, 'Number of people');
            hold(ax, 'on');
            
            if ~iscell(Plots)
                Plots = num2cell(Plots);
            end
            
            for pi = 1:length(Plots)
                plot(ax, XData, obj.getValue(Plots{pi}));
            end
            legend(ax, 'show');
            hold(ax, 'off');
        end
        
        function h = plotHeatMap(obj, varargin)
            Step = 1;
            
            for ii = 1:2:length(varargin)
                switch string(varargin{ii})
                    case "Step"
                        Step = varargin{ii + 1};
                end
            end
            
            h = obj.State(Step).plotHeatMap(varargin{:});
            h.Title = strcat(h.Title, " at ", datestr(obj.Date(Step), 'dd mmmm yyyy'));
        end
        
    end
    
end

