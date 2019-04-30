function [problem, state] = equationsThreePhaseBlackOilPolymer2(state0, state, model, dt, drivingForces, varargin)

%
%
% SYNOPSIS:
%   function [problem, state] = equationsThreePhaseBlackOilPolymer(state0, state, model, dt, drivingForces, varargin)
%
% DESCRIPTION: 
%   Assemble the linearized equations for a blackoil system,
%   computing both the residuals and the Jacobians. Returns the result as an
%   instance of the class LinearizedProblem which can be solved using instances of
%   LinearSolverAD.
%
%   A description of the modeling equations can be found in the directory
%   ad-eor/docs.
%
% PARAMETERS:
%   state0        - State at previous times-step
%   state         - State at current time-step
%   model         - Model instance
%   dt            - time-step
%   drivingForces - Driving forces (boundary conditions, wells, ...)
%   varargin      - optional parameters
%
% RETURNS:
%   problem - Instance of LinearizedProblem
%   state   - Updated state variable (fluxes, mobilities and more can be
%             stored, the wellSol structure is also updated in case of control switching)
%
% EXAMPLE:
%
% SEE ALSO: LinearizedProblem, LinearSolverAD, equationsOilWater, OilWaterPolymerModel
%

% Get linearized problem for oil/water/polymer system with black oil-style
% properties

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

    opt = struct('Verbose', mrstVerbose, ...
                 'reverseMode', false,...
                 'resOnly', false,...
                 'iteration', -1);

    opt = merge_options(opt, varargin{:});

    % Shorter names for some commonly used parts of the model and forces.
    s = model.operators;
    f = model.fluid;
    G = model.G;
    W = drivingForces.W;
    % Currently we do not support senario without wells.
    assert(isempty(drivingForces.bc) && isempty(drivingForces.src));

    % Properties at current timestep
    [p, sW, sG, rs, rv, c, cmax, wellSol] = model.getProps(state, ...
                                                      'pressure', 'water', 'gas', 'rs', 'rv', 'polymer', 'polymermax', 'wellsol');

    % Properties at previous timestep
    [p0, sW0, sG0, rs0, rv0, c0, cmax0, wellSol0] = model.getProps(state0, ...
                                                      'pressure', 'water', 'gas', 'rs', 'rv', 'polymer', 'polymermax', 'wellsol');

    [wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

    %Initialization of primary variables ----------------------------------
    st  = model.getCellStatusVO(state,  1-sW-sG,   sW,  sG);
    st0 = model.getCellStatusVO(state0, 1-sW0-sG0, sW0, sG0);

    if model.disgas || model.vapoil
        % X is either Rs, Rv or Sg, depending on each cell's saturation status
        x = st{1}.*rs + st{2}.*rv + st{3}.*sG;
        gvar = 'x';
    else
        x = sG;
        gvar = 'sG';
    end

    if ~opt.resOnly
        if ~opt.reverseMode
            % define primary varible x and initialize
            [p, sW, x, c, wellVars{:}] = ...
                initVariablesADI(p, sW, x, c, wellVars{:});
        else
            x0 = st0{1}.*rs0 + st0{2}.*rv0 + st0{3}.*sG0;
            % Set initial gradient to zero
            wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
            [p0, sW0, x0, c0, wellVars0{:}] = ...
                initVariablesADI(p0, sW0, x0, c0, wellVars0{:}); %#ok
            clear zw;
            [sG0, rs0, rv0] = calculateHydrocarbonsFromStatusBO(model, st0, 1-sW, x0, rs0, rv0, p0);
        end
    end

    if ~opt.reverseMode
        % Compute values from status flags. If we are in reverse mode, these
        % values have already converged in the forward simulation.
        [sG, rs, rv, rsSat, rvSat] = calculateHydrocarbonsFromStatusBO(model, st, 1-sW, x, rs, rv, p);
    end

    % We will solve for pressure, water and gas saturation (oil saturation
    % follows via the definition of saturations), polymer concentration and well rates + bhp.
    primaryVars = {'pressure', 'sW', gvar, 'polymer', wellVarNames{:}};

    % Evaluate relative permeability
    sO  = 1 - sW  - sG;
    sO0 = 1 - sW0 - sG0;
    if model.water
        sat = {sW, sO, sG};
        sat0 = {sW0, sO0, sG0};
    else
        sat = {sO, sG};
        sat0 = {sO0, sG0};
    end

    % Update state with AD-variables
    state = model.setProps(state  , {'s', 'pressure', 'rs', 'rv', 'polymer'}, {sat , p , rs , rv, c});
    state0 = model.setProps(state0, {'s', 'pressure', 'rs', 'rv'}, {sat0, p0, rs0, rv0});
    % Set up properties
    state = model.initPropertyContainers(state);
    pBH = wellVars{wellMap.isBHP};
    [b, pv] = model.getProps(state, 'ShrinkageFactors', 'PoreVolume');
    [b0, pv0] = model.getProps(state0, 'ShrinkageFactors', 'PoreVolume');
    [phaseFlux, flags] = model.getProps(state, 'PhaseFlux',  'PhaseUpwindFlag');

    [bW, bO, bG]       = deal(b{:});
    [bW0, bO0, bG0]    = deal(b0{:});
    [vW, vO, vG]       = deal(phaseFlux{:});
    [upcw, upco, upcg] = deal(flags{:});

    [pressures, mob, rho] = model.getProps(state, 'PhasePressures', 'Mobility', 'Density');
    [mobW, mobO, mobG] = deal(mob{:});
    bWvW = s.faceUpstr(upcw, bW).*vW;
    bOvO = s.faceUpstr(upco, bO).*vO;
    bGvG = s.faceUpstr(upcg, bG).*vG;
    % The first equation is the conservation of the water phase. This equation is
    % straightforward, as water is assumed to remain in the aqua phase in the
    % black oil model.
    water = (1/dt).*(pv.*bW.*sW - pv0.*bW0.*sW0);
    divWater = s.Div(bWvW);

    % Second equation: mass conservation equation for the oil phase at surface
    % conditions. This is any liquid oil at reservoir conditions, as well as
    % any oil dissolved into the gas phase (if the model has vapoil enabled).
    if model.vapoil
        % The model allows oil to vaporize into the gas phase. The conservation
        % equation for oil must then include the fraction present in the gas
        % phase.
        rvbGvG = s.faceUpstr(upcg, rv).*bGvG;
        % Final equation
        oil = (1/dt).*( pv .*(bO.* sO  + rv.* bG.* sG) - ...
                        pv0.*(bO0.*sO0 + rv0.*bG0.*sG0));
        divOil = s.Div(bOvO + rvbGvG);
    else
        oil = (1/dt).*(pv.*bO.*sO - pv0.*bO0.*sO0 );
        divOil = s.Div(bOvO);
    end

    % Conservation of mass for gas. Again, we have two cases depending on
    % whether the model allows us to dissolve the gas phase into the oil phase.
    if model.disgas
        % The gas transported in the oil phase.
        rsbOvO = s.faceUpstr(upco, rs).*bOvO;
        
        gas = (1/dt).*(pv.* (bG.* sG  + rs.* bO.* sO) - ...
                       pv0.*(bG0.*sG0 + rs0.*bO0.*sO0 ));
        divGas = s.Div(bGvG + rsbOvO);
    else
        gas = (1/dt).*(pv.*bG.*sG - pv0.*bG0.*sG0 );
        divGas = s.Div(bGvG);
    end

    % polymer in water equation :
    poro =  s.pv./G.cells.volumes;
    polymer = (s.pv.*(1-f.dps)/dt).*(pvMult.*bW.*sW.*c - ...
                                     pvMult0.*f.bW(p0).*sW0.*c0) + (s.pv/dt).* ...
              ( f.rhoR.*((1-poro)./poro).*(ads - ads0)) + s.Div(bWvP);

    % Applying correction to the polymer equation when the Jacobian is
    % prolematic for some cells.
    % Typically it is due to totally and almost non-existence of water.
    if ~opt.resOnly
        epsilon = 1.e-8;
        % the first way is based on the diagonal values of the resulting
        % Jacobian matrix
        eps = sqrt(epsilon)*mean(abs(diag(polymer.jac{4})));
        % sometimes there is no water in the whole domain
        if (eps == 0.)
            eps = epsilon;
        end
        % bad marks the cells prolematic in evaluating Jacobian
        bad = abs(diag(polymer.jac{4})) < eps;
        % the other way is to choose based on the water saturation
        polymer(bad) = c(bad);
    end

    eqs = {water, oil, gas, polymer};
    names = {'water', 'oil', 'gas', 'polymer'};
    types = {'cell', 'cell', 'cell', 'cell'};

    rho = {rhoW, rhoO, rhoG};
    mob = {mobW, mobO, mobG};
    sat = {sW, sO, sG};
    dissolved = model.getDissolutionMatrix(rs, rv);

    [eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                   {pW, p, pG}, sat, mob, rho, ...
                                                   dissolved, {c}, ...
                                                   drivingForces);

    % Finally, add in and setup well equations
    wc    = vertcat(W.cells);
    perf2well = getPerforationToWellMapping(W);
    sgn = vertcat(W.sign);

    wc_inj = wc(sgn(perf2well) > 0);
    cw     = c(wc_inj);

    % remove the old viscosity and applying the fully mixed viscosity
    muWMultW = muWMult(wc_inj);
    muWFullyMixed = model.fluid.muWMult(cw);

    mob{1}(wc_inj) = mob{1}(wc_inj) ./ muWFullyMixed .* muWMultW;


    if model.usingShear || model.usingShearLog || model.usingShearLogshrate
        % applying the shear effects
        mob{1}(wc) = mob{1}(wc)./shearMultW;
    end

    [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, p, mob, rho, dissolved, {c}, dt, opt);
    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

end

%--------------------------------------------------------------------------

function [wPoly, wciPoly, iInxW] = getWellPolymer(W)
    if isempty(W)
        wPoly = [];
        wciPoly = [];
        iInxW = [];
        return
    end
    inj   = vertcat(W.sign)==1;
    polInj = cellfun(@(x)~isempty(x), {W(inj).c});
    wPoly = zeros(nnz(inj), 1);
    W_inj = W(inj);
    wPoly(polInj) = vertcat(W_inj(polInj).c);
    wciPoly = rldecode(wPoly, cellfun(@numel, {W_inj.cells}));

    % Injection cells
    nPerf = cellfun(@numel, {W.cells})';
    nw    = numel(W);
    perf2well = rldecode((1:nw)', nPerf);
    compi = vertcat(W.compi);
    iInx  = rldecode(inj, nPerf);
    iInx  = find(iInx);
    iInxW = iInx(compi(perf2well(iInx),1)==1);
end

%--------------------------------------------------------------------------

% Effective adsorption, depending of desorption or not
function y = effads(c, cmax, f)
    if f.adsInx == 2
        y = f.ads(max(c, cmax));
    else
        y = f.ads(c);
    end
end

