function cfl = estimateSaturationCFL(model, state, dt, varargin)
    opt = struct('forces', []);
    opt = merge_options(opt, varargin{:});
    if isfield(state, 'flux')
        v = state.flux(model.operators.internalConn, :);
    else
        v = value(model.getProps(state, 'PhaseFlux'));
    end
    pv = model.operators.pv;
    
    [F, Q] = getFractionalFlowMagnitude(model, state);
    
    vT = sum(v, 2);
    
    xflow = ~(all(v >= 0, 2) | all(v <= 0, 2));
    
    l = model.operators.N(:, 1);
    r = model.operators.N(:, 2);
    
    flag = vT > 0;
    df_face = upstream(model, F, flag, xflow, l, r);
    
    rate_face = abs(vT).*abs(df_face);
    if ~isempty(Q)
        T = model.operators.T;
        cap_face = upstream(model, Q, flag, xflow, l, r);
        rate_face = rate_face + 2.*T.*cap_face;
    end
    
    nc = model.G.cells.num;
    % Accumulate into cell if flow is outgoing, or we have any kind of
    % cross-flow.
    rate_cell = accumarray(l, rate_face.*( flag | xflow), [nc, 1]) +...
                accumarray(r, rate_face.*(~flag | xflow), [nc, 1]);
    if ~isempty(opt.forces)
        if isfield(state.wellSol, 'flux') % Wells
            wflux = sum(vertcat(state.wellSol.flux), 2);
            wc = vertcat(opt.forces.W.cells);
            rate_cell(wc) = rate_cell(wc) + abs(wflux).*F(wc);
        end
    end
    cfl = (dt./pv).*rate_cell;
end

function df_face = upstream(model, F, flag, xflow, l, r)
    df_face = model.operators.faceUpstr(flag, F);
    df_face(xflow) = max(abs(F(l(xflow))), abs(F(r(xflow))));
end