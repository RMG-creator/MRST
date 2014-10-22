classdef ReservoirModel < PhysicalModel
%Base class for physical models
%
% SYNOPSIS:
%   model = ReservoirModel(G, rock, fluid)
%
% DESCRIPTION:
%   Extension of PhysicalModel class to accomodate reservoir-specific
%   features such as fluid and rock as well as commonly used phases and
%   variables.
%
% REQUIRED PARAMETERS:
%   G     - Simulation grid.
%
%   rock  - Valid rock used for the model.
%
%   fluid - Fluid model used for the model.
%
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   See class properties.
%
% RETURNS:
%   Class instance.
%
% SEE ALSO:
%   ThreePhaseBlackOilModel, TwoPhaseOilWaterModel, PhysicalModel

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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

    properties
        % The fluid model
        fluid
        % The rock (perm/poro/ntg)
        rock
        
        % Maximum relative pressure change
        dpMaxRel
        dpMaxAbs
        % Maximum absolute saturation change
        dsMaxRel
        dsMaxAbs
        % Water phase present
        water
        % Gas phase present
        gas
        % Oil phase present
        oil
        % Names of primary variables interpreted as saturations, i.e. so
        % that they will sum to one when updated.
        saturationVarNames
        % Names of well fields that may be updated by the model.
        wellVarNames
        
        % Use alternate tolerance scheme
        useCNVConvergence
        
        % CNV tolerance (inf-norm-like)
        toleranceCNV;
        
        % MB tolerance values (2-norm-like)
        toleranceMB;
        % Well tolerance if CNV is being used
        toleranceWellBHP;
        % Well tolerance if CNV is being used
        toleranceWellRate;
        
        % Input data used to instantiate the model
        inputdata
        % Add extra output to wellsol/states for derived quantities
        extraStateOutput
        % Output fluxes
        outputFluxes
    end
    
    methods
        function model = ReservoirModel(G, rock, fluid, varargin)
            model = model@PhysicalModel(G);
            
            model.dpMaxRel = inf;
            model.dpMaxAbs = inf;
            
            model.dsMaxAbs = .2;
            model.dsMaxRel = inf;
            
            model.nonlinearTolerance = 1e-6;
            model.inputdata = [];
            
            model.useCNVConvergence = false;
            model.toleranceCNV = 1e-3;
            model.toleranceMB = 1e-7;
            model.toleranceWellBHP = 1*barsa;
            model.toleranceWellRate = 1/day;
            
            model.saturationVarNames = {'sw', 'so', 'sg'};
            model.wellVarNames = {'qWs', 'qOs', 'qGs', 'bhp'};
            
            model.extraStateOutput = false;
            model.outputFluxes = true;
            
            model = merge_options(model, varargin{:});
            
            % Base class does not support any phases
            model.water = false;
            model.gas = false;
            model.oil = false;
            
            % Physical model
            model.fluid = fluid;
            model.rock  = rock;
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            % Generic update function for reservoir models containing wells

            % Split variables into three categories: Regular/rest variables, saturation
            % variables (which sum to 1 after updates) and well variables (which live
            % in wellSol and are in general more messy to work with).
            [restVars, satVars, wellVars] = model.splitPrimaryVariables(problem.primaryVariables);
            
            % Update saturations in one go
            state  = model.updateSaturations(state, dx, problem, satVars);
            
            if ~isempty(restVars)
                % Handle pressure seperately
                state = model.updateStateFromIncrement(state, dx, problem, 'pressure', model.dpMaxRel, model.dpMaxAbs);
                restVars = model.stripVars(restVars, 'pressure');

                % Update remaining variables (tracers, temperature etc)
                for i = 1:numel(restVars);
                     state = model.updateStateFromIncrement(state, dx, problem, restVars{i});
                end
            end

            % Update the wells
            if isfield(state, 'wellSol')
                state.wellSol = model.updateWellSol(state.wellSol, problem, dx, drivingForces, wellVars);
            end
            report = [];
        end

        function model = setupOperators(model, G, rock, varargin)
            % Set up divergence/gradient/transmissibility operators
            if isempty(G) || isempty(rock)
                warning('mrst:ReservoirModel', ...
                'Invalid grid/rock pair supplied. Operators have not been set up.')
                return;
            end
            model.operators = setupOperatorsTPFA(G, rock, varargin{:});
        end
        
        function [convergence, values] = checkConvergence(model, problem, varargin)
            if model.useCNVConvergence
                % Use convergence model similar to commercial simulator
                [conv_cells, v_cells] = CNV_MBConvergence(model, problem);
                [conv_wells, v_wells] = checkWellConvergence(model, problem);
                
                convergence = all(conv_cells) && all(conv_wells);
                values = [v_cells, v_wells];
            else
                % Use strict tolerances on the residual without any 
                % fingerspitzengefuhlen by calling the parent class
                [convergence, values] = checkConvergence@PhysicalModel(model, problem, varargin{:});
            end            
        end
        
        function [vararg, driving] = getDrivingForces(model, control) %#ok
            % Setup and pass on driving forces
            vararg = {};
            driving = struct('Wells', [], 'bc', [], 'src', []);
            
            if isfield(control, 'W') && ~isempty(control.W)
                vararg = [vararg, 'Wells', control.W];
                driving.Wells = control.W;
            end

            if isfield(control, 'bc') && ~isempty(control.bc)
                vararg = [vararg, 'bc', control.bc];
                driving.bc = control.bc;
            end
            
            if isfield(control, 'src') && ~isempty(control.src)
                vararg = [vararg, 'src', control.src];
                driving.src = control.src;
            end
        end
        
        function [fn, index] = getVariableField(model, name)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            switch(lower(name))
                case {'t', 'temperature'}
                    fn = 'T';
                    index = 1;
                case {'sw', 'water'}
                    index = find(strcmpi(model.saturationVarNames, 'sw'));
                    fn = 's';
                case {'so', 'oil'}
                    index = find(strcmpi(model.saturationVarNames, 'so'));
                    fn = 's';
                case {'sg', 'gas'}
                    index = find(strcmpi(model.saturationVarNames, 'sg'));
                    fn = 's';
                case {'s', 'sat', 'saturation'}
                    index = 1:numel(model.saturationVarNames);
                    fn = 's';
                case {'pressure', 'p'}
                    index = 1;
                    fn = 'pressure';
                case 'wellsol'
                    % Use colon to get all variables, since the wellsol may
                    % be empty
                    index = ':';
                    fn = 'wellSol';
                otherwise
                    % This will throw an error for us
                    [fn, index] = getVariableField@PhysicalModel(model, name);
            end
        end
        
        function [restVars, satVars, wellVars] = splitPrimaryVariables(model, vars)
            isSat   = cellfun(@(x) any(strcmpi(model.saturationVarNames, x)), vars);
            isWells = cellfun(@(x) any(strcmpi(model.wellVarNames, x)), vars);
            
            wellVars = vars(isWells);
            satVars  = vars(isSat);
                
            restVars = vars(~isSat & ~isWells);
        end
        
        function [isActive, phInd] = getActivePhases(model)
            isActive = [model.water, model.oil, model.gas];
            if nargout > 1
                phInd = find(isActive);
            end
        end
        
        function phNames = getPhaseNames(model)
            tmp = 'WOG';
            active = model.getActivePhases();
            phNames = tmp(active);
        end
        
        function i = getPhaseIndex(model, phasename)
            active = model.getPhaseNames();
            i = find(active == phasename);
        end 
        
        function state = updateSaturations(model, state, dx, problem, satVars)
            % Update saturations (likely state.s) under the constraint that
            % the sum of volume fractions is always equal to 1. This
            % assumes that we have solved for n - 1 phases when n phases
            % are present.
            if nargin < 5
                % Get the saturation names directly from the problem
                [~, satVars] = ...
                    splitPrimaryVariables(model, problem.primaryVariables);
            end
            if isempty(satVars)
                % No saturations passed, nothing to do here.
                return
            end
            % Solution variables should be saturations directly, find the missing
            % link
            saturations = lower(model.saturationVarNames);
            fillsat = setdiff(saturations, lower(satVars));
            assert(numel(fillsat) == 1)
            fillsat = fillsat{1};

            % Fill component is whichever saturation is assumed to fill up the rest of
            % the pores. This is done by setting that increment equal to the
            % negation of all others so that sum(s) == 0 at end of update
            solvedFor = ~strcmpi(saturations, fillsat);
            ds = zeros(model.G.cells.num, numel(saturations));

            tmp = 0;
            for i = 1:numel(saturations)
                if solvedFor(i)
                    v = model.getIncrement(dx, problem, saturations{i});
                    ds(:, i) = v;
                    % Saturations added for active variables must be subtracted
                    % from the last phase
                    tmp = tmp - v;
                end
            end
            ds(:, ~solvedFor) = tmp;
            % We update all saturations simultanously, since this does not bias the
            % increment towards one phase in particular.
            state = model.updateStateFromIncrement(state, ds, problem, 's', model.dsMaxRel, model.dsMaxAbs);
            
            % Ensure that values are within zero->one interval, and
            % re-normalize if any values were capped
            bad = any((state.s > 1) | (state.s < 0), 2);
            if any(bad)
                state.s(bad, :) = min(state.s(bad, :), 1);
                state.s(bad, :) = max(state.s(bad, :), 0);
                state.s(bad, :) = bsxfun(@rdivide, state.s(bad, :), sum(state.s(bad, :), 2));
            end
        end
        
        function wellSol = updateWellSol(model, wellSol, problem, dx, drivingForces, wellVars) %#ok
            % Update the wellSol struct
            if numel(wellSol) == 0
                % Nothing to be done if there are no wells
                return
            end
            if nargin < 6
                % Get the well variables directly from the problem,
                % otherwise assume that they are known by the user
                [~, ~, wellVars] = ...
                    splitPrimaryVariables(model, problem.primaryVariables);
            end
            
            for i = 1:numel(wellVars)
                wf = wellVars{i};
                dv = model.getIncrement(dx, problem, wf);

                if strcmpi(wf, 'bhp')
                    % Bottom hole is a bit special - we apply the pressure update
                    % limits here as well.
                    bhp = vertcat(wellSol.bhp);
                    dv = model.limitUpdateRelative(dv, bhp, model.dpMaxRel);
                    dv = model.limitUpdateAbsolute(dv, model.dpMaxAbs);
                end

                for j = 1:numel(wellSol)
                    wellSol(j).(wf) = wellSol(j).(wf) + dv(j);
                end
            end
        end
        
        function state = setPhaseData(model, state, data, fld, subs)
            % Given a structure field name and a cell array of data for
            % each phase, store them as columns with the given field name
            if nargin == 4
                subs = ':';
            end
            isActive = model.getActivePhases();
            
            ind = 1;
            for i = 1:numel(data)
                if isActive(i)
                    state.(fld)(subs, ind) = data{i};
                    ind = ind + 1;
                end
            end
        end
        
        function state = storeFluxes(model, state, vW, vO, vG)
            isActive = model.getActivePhases();
            
            internal = model.operators.internalConn;
            state.flux = zeros(numel(internal), sum(isActive));
            phasefluxes = {double(vW), double(vO), double(vG)};
            state = model.setPhaseData(state, phasefluxes, 'flux', internal);
        end
        
        function state = storeMobilities(model, state, mobW, mobO, mobG)
            isActive = model.getActivePhases();
            
            state.mob = zeros(model.G.cells.num, sum(isActive));
            mob = {double(mobW), double(mobO), double(mobG)};
            state = model.setPhaseData(state, mob, 'mob');
        end
        
        function state = storeUpstreamIndices(model, state, upcw, upco, upcg)
            isActive = model.getActivePhases();
            
            nInterfaces = size(model.operators.N, 1);
            state.upstreamFlag = false(nInterfaces, sum(isActive));
            mob = {upcw, upco, upcg};
            state = model.setPhaseData(state, mob, 'upstreamFlag');
        end
        
        function state = storebfactors(model, state, bW, bO, bG)
            isActive = model.getActivePhases();
            
            state.mob = zeros(model.G.cells.num, sum(isActive));
            b = {double(bW), double(bO), double(bG)};
            state = model.setPhaseData(state, b, 'bfactor');
        end
    end

    methods (Static)

    end

end

