function [x, y, K, Z_L, Z_V, L, values] = newtonFugacityEquilibrium(model, P, T, z, K, L)
% Single step of the newton update based on fugacity
    
    ncomp = model.fluid.getNumberOfComponents();
    ncell = numel(L);

    x = model.computeLiquid(L, K, z);
    y = model.computeVapor(L, K, z);
    
    [xAD, yAD] = deal(cell(1, ncomp));
    [xAD{:}, yAD{:}, L_ad] = initVariablesADI(x{:}, y{:}, L);
    [eqs, f_L, f_V, Z_L, Z_V] = model.equationsEquilibrium(P, T, xAD, yAD, z, L_ad, [], []);

    eq = cat(eqs{:});
    
    J = eq.jac{1};
    r = eq.val;
    upd = -J\r;
    cap = @(x) min(max(x, 0), 1);

    [dx, dy] = deal(cell(1, ncomp));
    
    dL = upd(end-ncell + 1: end);
    bad = isinf(dL);
    dL(bad) = sign(dL(bad));
    
    w = getRelax(L, dL, cap);
    for i = 1:ncomp
        subs = (ncell*(i-1) + 1):(ncell*i);
        dx{i} = upd(subs);
        dy{i} = upd(subs + ncomp*ncell);
        
        dx{i}(L == 0) = 0;
        dy{i}(L == 1) = 0;
        
        w_x = getRelax(x{i}, dx{i}, cap);
        w_y = getRelax(y{i}, dy{i}, cap);
        w = min(w, min(w_x, w_y));        
    end
    
    f_r = cell(1, ncomp);
    L = L + w.*dL;

    pure = L == 1 | L == 0;

    K = cell(1, ncomp);
    for i = 1:ncomp
        zmissing = z{i} == 0;
        x{i} = x{i} + w.*dx{i};
        y{i} = y{i} + w.*dy{i};
        x{i}(pure) = z{i}(pure);
        y{i}(pure) = z{i}(pure);
        
        f_r{i} = double(f_L{i})./double(f_V{i});
        f_r{i}(pure | zmissing) = 1;
        
        K{i} = y{i}./x{i};
        
        xmissing = x{i} == 0;
        gone = zmissing | xmissing | pure;
        K{i}(gone) = 1;
    end
    fn = @(v) cellfun(@(x) all(isfinite(x)), v);
    assert(all(fn(K)));
    assert(all(fn(x)));
    assert(all(fn(y)));
    assert(all(isfinite(L)));
    
    values = abs([f_r{:}] - 1);
end


function w = getRelax(x, dx, cap)
    w = (cap(x + dx) - x)./(dx);
    w(isnan(w)) = 0;
end