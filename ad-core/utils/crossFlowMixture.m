function compi = crossFlowMixture(flux, compi, map, conserveMass, net_injection_given)
    if nargin < 4
        conserveMass = false;
    end
    % Into wellbore
    flux_in = -min(flux, 0);
    if all(flux_in == 0)
        return
    end

    % Find net flux - this is what is possibily injected from the
    % surface connection with the top composition
    total_flux = sum(flux, 2);
    into_well_from_reservoir = sum_perf(-min(flux, 0), map);
    from_well_into_reservoir = sum_perf(max(flux, 0), map);
    net_flux = sum(from_well_into_reservoir, 2) - sum(into_well_from_reservoir, 2);
    
    net_injection = max(net_flux, 0);
    % Flux into well-bore plus net injection weighted with topside
    % composition
    top_in =  bsxfun(@times, net_injection, compi);
    if nargin > 4
        replace = ~isnan(net_injection_given);
        top_in(replace, :) = net_injection_given(replace).*compi(replace, :);
    end
    comp = into_well_from_reservoir + top_in;
    compT = sum(comp, 2);
    % Normalize to get fractions
    comp = bsxfun(@rdivide, comp, compT);
    if conserveMass
        % Ensure exact re-injection with "fractions" which are not in unit
        % range for injectors
        act = sum(top_in, 2) > 0;
        % Top net injected composition + sum of composition into well-bore,
        % divided by total flux out from well-bore
        out = sum(max(from_well_into_reservoir(act, :), 0), 2);
        tmp = (top_in(act, :) + into_well_from_reservoir(act, :))./max(out, 1e-12);
        comp(act, :) = tmp;
    end
    active = compT > 0;
    compi(active, :) = comp(active, :);
end

function out = sum_perf(v, map)
    nph = size(v, 2);
    nw = numel(map.W);
    out = zeros(nw, nph);
    for i = 1:nph
        out(:, i) = accumarray(map.perf2well, v(:, i));
    end
end