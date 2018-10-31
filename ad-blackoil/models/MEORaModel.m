classdef MEORaModel < TwoPhaseOilWaterModel
    % Oil/water/microbial system
    % This model is a two phase oil/water model, extended with the 
    % microbial phase in addition.
    % Microbe effects currently available for simulation:
    % Nothing
    
    properties
        % Substances in reservoir
        microbe
        nutrient
        metabolite
        biofilm
        stateplots
        % Variables
        yield_microbe
        yield_metabolite
        growth_max_microbe
        growth_max_metabolite
        halfsat_microbe
        halfsat_metabolite
        crit_val
        biosurf
        biopoly
        langmuir
    end
    
    methods
        function model = MEORaModel(G, rock, fluid, varargin)
            
            model = model@TwoPhaseOilWaterModel(G,rock, fluid);
            
            model.microbe = true;
            model.metabolite = true;
            model.nutrient = true;
            % pre-defining defaults
            model.yield_microbe = .5;
            model.yield_metabolite = .5;
            model.growth_max_microbe = .2/day;
            model.growth_max_metabolite = .2/day;
            model.halfsat_microbe = .5;
            model.halfsat_metabolite = .5;
            model.crit_val = 0;
            model.langmuir = [0, 0];
            
            % type/types of metabolite to be defined manually
            model.biosurf = false;
            model.biopoly = false;
            model.biofilm = false;
            modle.stateplots = false;
            % not available for more than 3 phases
            % oil/water/microbe/nutrient/metabolite
            model.outputFluxes = true;
            
            model.wellVarNames = {'qWs', 'qOs', 'qWMEOR', 'bhp'};
            
            model = merge_options(model, varargin{:});
        end
        
        
        function [problem, state] = getEquations(model, state0, state,...
                dt, drivingForces, varargin)
            if model.biofilm
            [problem, state] = equationsMEORbiofilm(state0, state, model,...
                dt, drivingForces, varargin{:});
            else
            [problem, state] = equationsMEORa(state0, state, model,...
                dt, drivingForces, varargin{:});
            end
        end
        
        function [fn, index] = getVariableField(model, name)
            % add metabolites here
            switch(lower(name))
                case 'microbe'
                    fn = 'm';
                    index = 1;
                case 'nutrient'
                    fn = 'n';
                    index = 1;
                case 'metabolite'
                    fn = 'meta';
                    index = 1;
                case 'biofilm'
                    fn = 'bio';
                    index = 1;
                otherwise
                    [fn, index] = getVariableField@TwoPhaseOilWaterModel(model, name);
            end
        end
        
        function [state, report] = updateState(model, state, problem, ...
                dx, drivingForces)
            [state, report] = updateState@TwoPhaseOilWaterModel(model, ...
               state, problem,  dx, drivingForces);
       
           if model.stateplots && rem(problem.iterationNo,5)==1
               xvals = linspace(0,1,model.G.cells.num);
       set(0, 'currentfigure', 1);
       plot(xvals, state.pressure);
       title(sprintf('Pressure'));
       
       set(0, 'currentfigure', 2);
       plot(xvals, state.s);
       title(sprintf('Saturation'));

       set(0, 'currentfigure', 3);
       plot(xvals, state.m);
       title(sprintf('Microbe concentration'));

       set(0, 'currentfigure', 4);
       plot(xvals, state.n);
       title(sprintf('Nutrient concentration'));
       
       set(0, 'currentfigure', 5);
       plot(xvals, state.meta);
       title(sprintf('Metabolite concentration'));
       if model.biofilm
        set(0, 'currentfigure', 6);
       plot(xvals, state.bio);
       title(sprintf('biofilm concentration'));
       end
       drawnow;
           end
        end
%         
%         function varargout = evaluteRelPerm(model, sat, varargin)
%             active = model.getActivePhases();
%             nph = sum(active);
%             assert(nph == numel(sat), ...
%             'The number of saturations must equal the number of active phases.')
%             varargout = cell(1, nph);
%             [varargout{:}] = model.fluid.relperm([sat{1} sat{2}]);
%         end
    end
end