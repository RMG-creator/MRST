classdef ThreePhaseSurfactantPolymerModel < ThreePhaseBlackOilModel
    %Three-phase black-oil model with support for surfactant and polymer injection
    %
    % SYNOPSIS:
    %   model = ThreePhaseSurfactantPolymerModel(G, rock, fluid, varargin)
    %
    % DESCRIPTION: 
    %   Fully implicit three phase blackoil model with polymer.
    %
    % PARAMETERS:
    %   G        - Grid
    %   rock     - Rock structure
    %   fluid    - Fluid structure
    %   varargin - optional parameters
    %
    % RETURNS:
    %   class instance
    %
    % EXAMPLE:
    %
    % SEE ALSO:  equationsThreePhaseBlackOilPolymer, OilWaterPolymerModel
    %

    properties
        % Surfactant and Polymer present
        polymer
        surfactant
    end

    methods
        function model = ThreePhaseSurfactantPolymerModel(G, rock, fluid, varargin)
            model = model@ThreePhaseBlackOilModel(G, rock, fluid, varargin{:});

            % This is the model parameters for oil/water/gas/surfactant/polymer system
            model.polymer = true;
            model.surfactant = true;
%           Add veloc calculation
            model = model.setupOperators(G, rock, varargin{:});
            model = merge_options(model, varargin{:});
        end

        % --------------------------------------------------------------------%
        function model = setupOperators(model, G, rock, varargin)
            model.operators.veloc = computeVelocTPFA(G, model.operators.internalConn);
            model.operators.sqVeloc = computeSqVelocTPFA(G, model.operators.internalConn);
        end

        % --------------------------------------------------------------------%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsThreePhaseSurfactantPolymer(state0, state, ...
                model, dt, drivingForces, varargin{:});
        end

        % --------------------------------------------------------------------%
        function state = validateState(model, state)
            state = validateState@ThreePhaseBlackOilPolymerModel(model, state);
            nc = model.G.cells.num;
            model.checkProperty(state, 'Polymer', [nc, 1], [1, 2]);
            model.checkProperty(state, 'Polymermax', [nc, 1], [1, 2]);
            model.checkProperty(state, 'Surfactant', [nc, 1], [1, 2]);
            model.checkProperty(state, 'SurfactantMax', [nc, 1], [1, 2]);
        end

        % --------------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            [state, report] = updateState@ThreePhaseBlackOilModel(model, ...
               state, problem,  dx, drivingForces);

            % c denotes concentration of polymer, cs is concentration of surfactant.    
            if model.polymer
                cp = model.getProp(state, 'polymer');
                cp = min(cp, model.fluid.cpmax);
                state = model.setProp(state, 'polymer', max(cp, 0));
            end
            if model.surfactant
                cs = model.getProp(state, 'surfactant');
                cs = min(cs, model.fluid.csmax);
                state = model.setProp(state, 'surfactant', max(cs, 0) );
            end
        end

        % --------------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@ThreePhaseBlackOilPolymerModel(model, state0, state, dt, drivingForces);

            if model.polymer
                cp     = model.getProp(state, 'polymer');
                cpmax  = model.getProp(state, 'polymermax');
                state = model.setProp(state, 'polymermax', max(cpmax, cp));
            end
            if model.surfactant
                cs     = model.getProp(state, 'surfactant');
                csmax  = model.getProp(state, 'surfactantmax');
                state = model.setProp(state, 'surfactantmax', max(csmax, cs));
            end
        end
        
        % --------------------------------------------------------------------%
        function model = validateModel(model, varargin)
            if isempty(model.FlowPropertyFunctions)
                model.FlowPropertyFunctions = SurfactantPolymerFlowPropertyFunctions(model);
            end
            if isempty(model.FluxDiscretization)
                model.FluxDiscretization = PolymerFluxDiscretization(model);
            end
            model = validateModel@ThreePhaseBlackOilModel(model, varargin{:});
        end

        % --------------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            switch(lower(name))
                case {'polymer', 'polymermax'}
                    cp = model.getComponentNames();
                    index = find(strcmpi(cp, 'polymer'));
                    if strcmpi(name, 'polymer')
                        fn = 'cp';
                    else
                        fn = 'cpmax';
                    end
                case {'surfactant', 'surfactantmax'}
                    cs = model.getComponentNames();
                    index = find(strcmpi(cs, 'surfactant'));
                    if strcmpi(name, 'surfactant')
                        fn = 'cs';
                    else
                        fn = 'csmax';
                    end    
                case 'qwpoly'
                    fn = 'qWPoly';
                    index = ':';
                case 'qwsft'
                    fn = 'qWSft';
                    index = ':';
                otherwise
                    [fn, index] = getVariableField@ThreePhaseBlackOilModel(...
                                    model, name, varargin{:});
            end
        end

        % --------------------------------------------------------------------%
        function names = getComponentNames(model)
            names = getComponentNames@ThreePhaseBlackOilModel(model);
            if model.polymer && model.surfactant
                names{end+1} = 'polymer';
                names{end+2} = 'surfactant';
            elseif model.polymer && ~model.surfactant
                names{end+1} = 'polymer';
            elseif model.surfactant
                names{end+1} = 'surfactant';
            end
        end
        
        % --------------------------------------------------------------------%
        function state = storeSurfData(model, state, s, cs, Nc, sigma)
            state.SWAT    = double(s);
            state.SURFACT = double(cs);
            state.SURFCNM = log(double(Nc))/log(10);
            state.SURFST  = double(sigma);
            % state.SURFADS = double(ads);
        end
        
        % --------------------------------------------------------------------%
        function scaling = getScalingFactorsCPR(model, problem, names, solver)
            nNames = numel(names);

            scaling = cell(nNames, 1);
            handled = false(nNames, 1);

            for iter = 1:nNames
                name = lower(names{iter});
                switch name
                    case 'polymer'
                        s = 0;
                    otherwise
                        continue
                end
                sub = strcmpi(problem.equationNames, name);

                scaling{iter} = s;
                handled(sub) = true;
            end
            if ~all(handled)
                % Get rest of scaling factors
                other = getScalingFactorsCPR@ThreePhaseBlackOilModel(model, problem, names(~handled), solver);
                [scaling{~handled}] = other{:};
            end
        end

        %% TODO CHECK WHAT IS THIS AND HOW TO CHANGE IT
        
        function [eq, src] = addComponentContributions(model, cname, eq, component, src, force)
        % For a given component conservation equation, compute and add in
        % source terms for a specific source/bc where the fluxes have
        % already been computed.
        %
        % PARAMETERS:
        %
        %   model  - (Base class, automatic)
        %
        %   cname  - Name of the component. Must be a property known to the
        %            model itself through `getProp` and `getVariableField`.
        %
        %   eq     - Equation where the source terms are to be added. Should
        %            be one value per cell in the simulation grid (model.G)
        %            so that the src.sourceCells is meaningful.
        %
        %   component - Cell-wise values of the component in question. Used
        %               for outflow source terms only.
        %
        %   src    - Source struct containing fields for fluxes etc. Should
        %            be constructed from force and the current reservoir
        %            state by `computeSourcesAndBoundaryConditionsAD`.
        %
        %   force  - Force struct used to produce src. Should contain the
        %            field defining the component in question, so that the
        %            inflow of the component through the boundary condition
        %            or source terms can accurately by estimated.
            if isempty(force)
                return
            end
            c = model.getProp(force, cname);
            cells = src.sourceCells;
            switch lower(cname)
              case {'polymer'}
                % Water based EOR, multiply by water flux divided by
                % density and add into corresponding equation
                qW = src.phaseMass{1}./model.fluid.rhoWS;
                isInj = qW > 0;
                qC = (isInj.*c + ~isInj.*component(cells)).*qW;
              otherwise
                error(['Unknown component ''', cname, '''. BC not implemented.']);
            end
            eq(cells) = eq(cells) - qC;
            src.components{end+1} = qC;
        end
        
        function [names, types] = getExtraWellEquationNames(model)
            [names, types] = getExtraWellEquationNames@ThreePhaseBlackOilModel(model);
            if model.polymer && model.surfactant
                names{end+1} = 'polymerWells';
                names{end+2} = 'surfactantWells';
                types{end+1} = 'perf';
                types{end+2} = 'perf';
            end
        end

        function names = getExtraWellPrimaryVariableNames(model)
            names = getExtraWellPrimaryVariableNames@ThreePhaseBlackOilModel(model);
            if model.polymer && model.surfactant
                names{end+1} = 'qWPoly';
                names{end+2} = 'qWSft';
            end
        end

        function [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration)
            [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions@ThreePhaseBlackOilModel(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration);
            if model.polymer && model.surfactant
                assert(model.water, 'Surfactant-Polymer injection requires a water phase.');
                f = model.fluid;

                % Water is always first
                wix = 1;
                cqWs = qMass{wix}./f.rhoWS; % connection volume flux at surface condition

                if well.isInjector()
                    concWellp = model.getProp(well.W, 'polymer');
                    concWells = model.getProp(well.W, 'surfactant');
                    cqP = concWellp.*cqWs;
                else
                    pixp = strcmpi(model.getComponentNames(), 'polymer');
                    pixs = strcmpi(model.getComponentNames(), 'surfactant');
                    concWellp = packed.components{pixp};
                    concWells = packed.components{pixs};

                    a = f.muWMult(f.cmax).^(1-f.mixPar);
                    cbarw     = concWellp/f.cmax;

                    % the term (a + (1 - a).*cbarw) account for the
                    % todd-longstaff mixing factor, which model the fact that for
                    % not-fully mixed polymer solution the polymer does not
                    % travel at the same velocity as water. See the governing
                    % equation for polymer (e.g. equationsOilWaterPolymer.m)
                    cqP = concWellp.*cqWs./(a + (1-a).*cbarw);
                end
                cqS = concWells.*cqWs;

                qwpoly = packed.extravars{strcmpi(packed.extravars_names, 'qwpoly')};
                qwsft = packed.extravars{strcmpi(packed.extravars_names, 'qwsft')};

                compEqs{end+1} = qwpoly - sum(concWellp.*cqWs);
                compEqs{end+2} = qwsft - sum(cqWs);
                compSrc{end+1} = cqP;
                compSrc{end+2} = cqS;
                eqNames{end+1} = 'polymerWells';
                eqNames{end+2} = 'surfactantWells';       
            end
        end
    end
end

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}