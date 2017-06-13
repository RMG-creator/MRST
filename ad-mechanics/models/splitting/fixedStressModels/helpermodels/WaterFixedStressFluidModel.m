classdef WaterFixedStressFluidModel < WaterModel
%
%
% SYNOPSIS:
%   model = WaterFixedStressFluidModel(G, rock, fluid, varargin)
%
% DESCRIPTION: Single phase fluid model to be used with fixed stress splitting
% method. The model handles the fluid equations of the splitting scheme
%
% PARAMETERS:
%   G        - Grid
%   rock     - rock structure
%   fluid    - fluid structure
%   varargin -
%
% RETURNS:
%   class instance
%
% EXAMPLE:
%
% SEE ALSO: WaterModel, MechFluidFixedStressSplitModel, MechFluidSplitModel
%
%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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

    methods
        function model = WaterFixedStressFluidModel(G, rock, fluid, varargin)
            model = model@WaterModel(G, rock, fluid);
            model = merge_options(model, varargin{:});
        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)
            % Setup the equations for the fluid. The drivingForce contains
            % the volumetric changes computed from the mechanical equations.

            opt = struct('Verbose', mrstVerbose, ...
                         'reverseMode', false,...
                         'resOnly', false,...
                         'iteration', -1);  % Compatibility only

            opt = merge_options(opt, varargin{:});

            [p, wellSol] = model.getProps(state, 'pressure', 'wellsol');

            [wellVars, wellVarNames, wellMap] = ...
                model.FacilityModel.getAllPrimaryVariables(wellSol);

            if ~opt.resOnly,
                [p, wellVars{:}] = initVariablesADI(p, wellVars{:});
            end

            fnew = drivingForces.fixedStressTerms.new;
            mechTerm.new = fnew.pTerm.*p - fnew.sTerm;

            fold = drivingForces.fixedStressTerms.old;
            mechTerm.old = fold.pTerm.*p - fold.sTerm;

            otherDrivingForces = rmfield(drivingForces, 'fixedStressTerms');

            [eqs, names, types, state] = equationsWaterMech(state0, p, ...
                                                            wellVars, state, ...
                                                            model, dt, ...
                                                            mechTerm, ...
                                                            otherDrivingForces, ...
                                                            'iteration', ...
                                                            opt.iteration);

            primaryVars = {'pressure', wellVarNames{:}};

            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end

        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@WaterModel(model);
            % divergence term
            % struct mechTerm.new and mechTerm.old
            forces.fixedStressTerms = [];
        end

        function fds = getAllVarsNames(model)
        % list of all the variable names that are used by this fluid model.
        fds = {'wellSol', 'pressure'};
        end

    end

end
