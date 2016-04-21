function checkComponentMassBalance(model, state0, states, schedule, n)
% Check mass balance of a simulator run and print to screen
    if nargin < 5
        n = numel(states);
    end
    
    states = states(1:n);
    dt = schedule.step.val(1:n);

    ws = cellfun(@(x) x.wellSol, states, 'UniformOutput', false);    
    ncomp = numel(states{1}.components);
    
    fracx = getMassFraction(states{end}.x, model.EOSModel.fluid);
    fracy = getMassFraction(states{end}.y, model.EOSModel.fluid);

    fracx0 = getMassFraction(state0.x, model.EOSModel.fluid);
    fracy0 = getMassFraction(state0.y, model.EOSModel.fluid);
    
    rhoO0 = model.EOSModel.computeDensity(state0.pressure, state0.x, state0.Z_L, state0.T);
    rhoG0 = model.EOSModel.computeDensity(state0.pressure, state0.y, state0.Z_V, state0.T);

    watOffset = model.water;
    
    for i = 1:ncomp
        wcomp = bsxfun(@times, getWellComponent(ws, i), dt);
        wcomp = wcomp(:);
        
        injected = sum(wcomp(wcomp > 0));
        produced = abs(sum(wcomp(wcomp < 0)));
        
        oil = sum(fracx(:, i).*states{end}.s(:, 1 + watOffset).*states{end}.rho(:, 1 + watOffset).*model.operators.pv);
        gas = sum(fracy(:, i).*states{end}.s(:, 2 + watOffset).*states{end}.rho(:, 2 + watOffset).*model.operators.pv);
        
        oil0 = sum(fracx0(:, i).*state0.s(:, 1 + watOffset).*rhoO0.*model.operators.pv);
        gas0 = sum(fracy0(:, i).*state0.s(:, 2 + watOffset).*rhoG0.*model.operators.pv);
        
        doil = oil - oil0;
        dgas = gas - gas0;
        
        fprintf('Component ''%s'':\n', model.EOSModel.fluid.names{i});
        fprintf('* Oil mass: From %1.2e -> %1.2e (net %+1.2e)\n', oil0, oil, doil);
        fprintf('* Gas mass: From %1.2e -> %1.2e (net %+1.2e)\n', gas0, gas, dgas);
        fprintf('* Injected: %1.2e\n', injected);
        fprintf('* Produced: %1.2e\n', produced);

        in = injected + oil0 + gas0;
        out = produced + oil + gas;
        fprintf('* Start: %1.2e, End %1.2e - %1.2f%%.\n', in, out, 100*in./out);
    end
end


function frac = getMassFraction(components, fluid)
    ncomp = numel(components);
    nv = numel(components{1});
    mass = zeros(nv, ncomp);
    for i = 1:numel(components)
        mass(:, i) = fluid.molarMass(i).*components{i};
    end
    frac = bsxfun(@rdivide, mass, sum(mass, 2));
end

function d = getWellComponent(ws, compix)
    nw = numel(ws{1});
    d = zeros(numel(ws), nw);
    for i = 1:nw
        d(:, i) = cellfun(@(w) w(i).components(compix), ws);
    end
end